# =============================================================================
# Simulation.jl
#
# Módulo Simulation: encadena los módulos de SIRIUS para simular un escenario
# completo. Es lo que ejecutan scripts/main.jl, el lanzador de VS Code,
# benchmark.jl y MonteCarlo.jl.
#
# Flujo de run_simulation (JSON completo -> resultados):
#   1. Configuración (Environment.load_config) y error de navegación de las cámaras
#   2. Propagación de órbitas, posiciones de las estaciones y dirección del Sol
#   3. Objetivos: forma (Geometry), actitud (Attitude) y magnitud vista desde cada
#      observador (Photometry)
#   4. Observadores: actitud, orientación de la cámara y simulación de la cámara (CameraSim)
#   5. Exportación de CSV (Environment)
#
# Funciones:
#   run_spec                   : spec -> JSON completo en config/generated/ -> simulación
#   run_simulation             : simulación a partir de un JSON completo
#   propagate_scenario         : paso 2
#   show_orbits                : figura 3D de las órbitas (si graphic_display_flag)
#   compute_targets_photometry : paso 3
#   simulate_observers!        : paso 4
#
# Dependencias: todos los módulos de src/, que Project.jl incluye antes que este.
# =============================================================================

module Simulation

export run_simulation, run_spec

using ..Types: AttitudeData, LOSMeasurement, TrajectoryData
using ..Constants: DATA_DIR
using ..Environment
using ..Propagator: propagate
using ..Attitude: bind_attitude, compute_attitude, camera_dcm_history
using ..Geometry: DebrisModel, generate_debris
using ..Photometry: compute_apparent_magnitude
using ..CameraSim: camera_sim!
using ..Visualization: setup_scene, add_orbit!
using ..Scenario: generate_config_file, read_json

using SPICE: utc2et
using StaticArrays

"""
    run_simulation(cfg_path, output_path) -> Nothing

Simula un escenario a partir de su JSON completo (el que escribe el generador en
`config/generated/`) y escribe los resultados en `output_path`.

Por qué existe: es el punto de entrada común de la simulación; `run_spec` la llama
tras generar el JSON y MonteCarlo.jl la llama directamente con JSON ya hechos.

# Argumentos
- `cfg_path::AbstractString`: JSON completo del escenario.
- `output_path::AbstractString`: carpeta de resultados (CSV y, si se piden, PNG).

# Devuelve
- `nothing`. Imprime el tiempo de ejecución.

# Excepciones
- Las de `load_config` (JSON no válido) y las de cada módulo de la simulación.
"""
function run_simulation(cfg_path::AbstractString, output_path::AbstractString)
    t0 = time()
    load_kernels(DATA_DIR)                              # efemérides SPICE (solo la primera vez)

    env = load_config(cfg_path)
    env.navigation_error_flag &&
        nav_err!(env.objects, env.obs_keys, env.obs_pos_sigma, env.obs_vel_sigma)

    data, sun_dir = propagate_scenario(env)
    env.graphic_display_flag && show_orbits(data)

    photometry_data = compute_targets_photometry(env, data, sun_dir)
    los_data = simulate_observers!(env, data, photometry_data, sun_dir, output_path)

    export_data_csvs!(data, los_data, env; output_root = output_path)
    println("Execution time: ", time() - t0, " s")
    return nothing
end

"""
    propagate_scenario(env) -> (Dict{Symbol,TrajectoryData}, Vector{SVector{3,Float64}})

Propaga todas las órbitas del escenario y calcula, en los mismos instantes
0:Δt:T, las posiciones de las estaciones y la dirección del Sol.

Por qué existe: reúne todo lo que depende solo del tiempo, antes de la actitud y
la fotometría.

# Argumentos
- `env`: configuración leída con `load_config`.

# Devuelve
- `(data, sun_dir)`: trayectorias por nombre (objetivos, observadores, `nav_*` y
  estaciones; ECI, km) y dirección unitaria del Sol en ECI en cada instante.

# Excepciones
- Las de `propagate` y SPICE (kernels no cargados).
"""
function propagate_scenario(env)
    duration_s = env.Duration_hours * 3600                      # duración (s)
    t_samples = 0.0:env.Data_time_step_seconds:duration_s       # misma malla que propagate

    data = propagate(env.objects, env.init_epoch, duration_s, env.Data_time_step_seconds)
    merge!(data, gs_positions(env.stations_llh, env.init_epoch, t_samples))
    sun_dir = sun_dir_eci.(Ref(utc2et(env.init_epoch)), t_samples)
    return data, sun_dir
end

"""
    show_orbits(data) -> Nothing

Muestra en una figura 3D (GLMakie) las órbitas propagadas, sin las trayectorias
con error de navegación (`nav_*`).

Por qué existe: visualización rápida del escenario cuando `graphic_display_flag` es `true`.

# Argumentos
- `data::Dict{Symbol,TrajectoryData}`: trayectorias de `propagate_scenario`.

# Devuelve
- `nothing`.

# Excepciones
- Ninguna propia (las de GLMakie si no hay pantalla).
"""
function show_orbits(data)
    fig, ax = setup_scene()
    for (name, traj) in data
        startswith(String(name), "nav_") && continue           # mismas órbitas con error de navegación
        add_orbit!(ax, traj.r, String(name))
    end
    display(fig)
    return nothing
end

"""
    compute_targets_photometry(env, data, sun_dir) -> Dict{Tuple{Symbol,Symbol},NamedTuple}

Para cada objetivo: genera su modelo fotométrico (bloque `shape`), calcula su
actitud (bloque `attitude`) y su magnitud aparente vista desde cada observador.

Por qué existe: agrupa todo lo que depende solo de los objetivos, antes de simular
las cámaras.

# Argumentos
- `env`: configuración leída con `load_config`.
- `data`, `sun_dir`: salida de `propagate_scenario`.

# Devuelve
- `Dict` (observador, objetivo) -> resultado de `compute_apparent_magnitude`.

# Excepciones
- Las de `generate_debris`, `compute_attitude` y `compute_apparent_magnitude`.
"""
function compute_targets_photometry(env, data, sun_dir)
    photometry_data = Dict{Tuple{Symbol,Symbol}, NamedTuple}()
    no_context = (camera_mount = nothing, objective_names = nothing, trajectories = nothing)

    for obj_name in env.obj_keys
        debris = generate_debris(env.obj_shape[obj_name])          # malla, áreas y ópticas
        spec = bind_attitude(env.obj_attitude[obj_name], no_context)   # los objetivos no tienen cámara
        attitude = compute_attitude(spec, data[obj_name], sun_dir)

        for obs_name in env.obs_keys
            photometry_data[(obs_name, obj_name)] = compute_apparent_magnitude(
                data[obs_name], data[obj_name], sun_dir, debris, env.photometry_config, attitude.DCM_b2i)
        end
    end
    return photometry_data
end

"""
    simulate_observers!(env, data, photometry_data, sun_dir, output_path) -> Dict

Para cada observador: calcula la actitud del satélite y la orientación de su
cámara (C_ci = C_cb · C_bi) y simula la cámara, que escribe imágenes y posiciones
en píxeles y devuelve las medidas LOS.

Por qué existe: agrupa todo lo que depende de los observadores; el modo `tracking`
necesita aquí las trayectorias de los objetivos.

# Argumentos
- `env`: configuración leída con `load_config`.
- `data`, `sun_dir`: salida de `propagate_scenario`.
- `photometry_data`: salida de `compute_targets_photometry`.
- `output_path::AbstractString`: carpeta de resultados.

# Devuelve
- `Dict{Symbol, Dict{Symbol, Vector{LOSMeasurement}}}`: medidas LOS por observador y objetivo.

# Excepciones
- Las de `bind_attitude` (p. ej. objetivo de `tracking` inexistente) y de `camera_sim!`.
"""
function simulate_observers!(env, data, photometry_data, sun_dir, output_path)
    los_data = Dict{Symbol, Dict{Symbol, Vector{LOSMeasurement}}}()

    # Trayectorias de los objetivos por su nombre del JSON (sin "obj_"), para el modo "tracking"
    obj_trajs = Dict(String(k)[length("obj_")+1:end] => data[k] for k in env.obj_keys)

    for obs_name in env.obs_keys
        C_cb = env.obs_mount[obs_name]                             # montaje de la cámara (cuerpo -> cámara)
        context = (camera_mount = C_cb, objective_names = keys(obj_trajs), trajectories = obj_trajs)
        spec = bind_attitude(env.obs_attitude[obs_name], context)
        C_ci = camera_dcm_history(C_cb, compute_attitude(spec, data[obs_name], sun_dir))   # ECI -> cámara
        camera_sim!(data, los_data, photometry_data, C_ci, obs_name, env.obj_keys, output_path, env.camera_config)
    end
    return los_data
end

"""
    run_spec(spec_path, output_path=nothing) -> Nothing

Ejecuta la simulación completa a partir de una spec: genera el JSON completo en
`config/generated/<nombre de la spec>.json` y simula con él.

Por qué existe: es lo que se usa normalmente (terminal, VS Code, REPL); la spec es
lo único que se escribe a mano.

# Argumentos
- `spec_path::AbstractString`: spec escrita a mano (normalmente en `config/specs/`).
- `output_path`: carpeta de resultados; con `nothing` (por defecto), `output/<nombre de la spec>/`.

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si el archivo no existe, no es una spec (p. ej. un JSON completo
  de `config/generated/`) o la spec no es válida.
- Cualquier error de la simulación.

# Ejemplo
```julia
run_spec("config/specs/ejemplo_actitud.json")
```
"""
function run_spec(spec_path::AbstractString, output_path::Union{AbstractString,Nothing} = nothing)
    spec_path = abspath(spec_path)
    name = splitext(basename(spec_path))[1]                # nombre de la spec, sin .json
    root = normpath(joinpath(@__DIR__, ".."))              # raíz del repositorio

    # Sin carpeta indicada: una por escenario, para no mezclar resultados de specs distintas
    output_path === nothing && (output_path = joinpath(root, "output", name))

    # Un JSON completo tiene los elementos (SMA...) en cada objetivo; una spec tiene "orbit"/"pattern"
    spec = read_json(spec_path)
    objs = get(spec, "objectives", Dict{String,Any}())
    any(o -> o isa AbstractDict && haskey(o, "SMA"), values(objs)) && error(
        "'$spec_path' is a full scenario JSON, not a spec. Pass the spec from config/specs/ " *
        "(or call run_simulation(path, output_path) to simulate a full JSON).")

    # 1. Spec -> JSON completo en config/generated/ (mismo nombre que la spec)
    cfg_path = joinpath(root, "config", "generated", name * ".json")
    generate_config_file(spec_path, cfg_path)

    # 2. Simulación con el JSON generado
    run_simulation(cfg_path, normpath(abspath(output_path)))
    return nothing
end

end # module Simulation
