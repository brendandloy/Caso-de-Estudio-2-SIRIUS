# =============================================================================
# main.jl
#
# Script principal de SIRIUS: ejecuta la simulación completa a partir de una spec.
#
# Uso:
#   julia --project=. scripts/main.jl <spec.json> [<output_path>]
#   (desde VS Code: scripts/run_simulation_vscode.jl, cambiando la variable `spec`)
#
# Pasos:
#   1. Genera el JSON completo del escenario en config/generated/<nombre de la spec>.json
#      (si ya existe, se sobrescribe; guarda las semillas usadas, así que se puede repetir).
#   2. Simula con ese JSON: propagación -> forma y actitud de los objetivos -> fotometría ->
#      actitud de los observadores y orientación de su cámara -> cámara -> CSV.
# Sin <output_path>, los resultados se escriben en output/<nombre de la spec>/.
#
# Funciones:
#   run_spec       : spec -> JSON completo -> simulación (punto de entrada)
#   run_simulation : simulación a partir de un JSON completo (la usa también MonteCarlo.jl)
# =============================================================================

push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Project
using Project.Types
using Project.Propagator
using Project.Attitude
using Project.Scenario       # generador: spec -> JSON completo
using Project.Environment
using Project.Visualization
using Project.CameraSim
using Project.Geometry
using Project.Photometry
using Project.Constants: DATA_DIR    # carpeta data/ de la raíz (kernels SPICE)

using Dates
using SPICE
using StaticArrays

load_kernels(DATA_DIR)

#--------------------------------------CONFIG-------------------------------------------------



function run_simulation(
    cfg_path::AbstractString,
    output_path::AbstractString
)
    t0 = time()

    env = load_config(cfg_path)

    if env.navigation_error_flag
    nav_err!(env.objects, env.obs_keys, env.obs_pos_sigma, env.obs_vel_sigma)
    end

    #------------------------------------MAIN SCRIPT----------------------------------------------

    Time_samples = 0.0:env.Data_time_step_seconds:env.Duration_hours*3600

    data = propagate(env.objects, env.init_epoch, env.Duration_hours*3600, env.Data_time_step_seconds) 
    gs_data = gs_positions(env.stations_llh, env.init_epoch, Time_samples)
    sun_dir = sun_dir_eci.(Ref(utc2et(env.init_epoch)), Time_samples)

    # Merge everything
    merge!(data, gs_data)

    # Visualización
    if env.graphic_display_flag
        fig, ax = setup_scene()
        for (name, traj) in data
            startswith(String(name), "nav_") && continue
            add_orbit!(ax, traj.r, String(name))
        end
        display(fig)
    end


    # Modelo fotométrico de cada objetivo (malla, áreas y ópticas) según su bloque "shape"
    debris_models = Dict{Symbol, Dict{Symbol,Any}}()
    for obj_name in env.obj_keys
        debris_models[obj_name] = generate_debris(env.obj_shape[obj_name])
    end

    attitude_data = Dict{Symbol,AttitudeData}()
    photometry_data = Dict{Tuple{Symbol,Symbol}, NamedTuple}()
    los_data = Dict{Symbol, Dict{Symbol, Vector{LOSMeasurement}}}()
        
    for obj_name in env.obj_keys
        # Actitud según el modo del bloque "attitude" (spin, random, fixed, pointing...); sin cámara
        spec = bind_attitude(env.obj_attitude[obj_name], (camera_mount = nothing, objective_names = nothing, trajectories = nothing))
        attitude_data[obj_name] = compute_attitude(spec, data[obj_name], sun_dir)
        for obs_name in env.obs_keys
            photometry_data[(obs_name, obj_name)] = compute_apparent_magnitude(
                data[obs_name], data[obj_name], sun_dir,  debris_models[obj_name], env.photometry_config, attitude_data[obj_name].DCM_b2i
            )
        end
    end

    # Trayectorias de los objetivos por su nombre del JSON (sin "obj_"), para el modo "tracking"
    obj_trajs = Dict(String(k)[length("obj_")+1:end] => data[k] for k in env.obj_keys)

    for obs_name in env.obs_keys
        # Actitud del satélite observador; bind_attitude da a "tracking" la cámara y los objetivos
        C_cb = env.obs_mount[obs_name]
        spec = bind_attitude(env.obs_attitude[obs_name],
                             (camera_mount = C_cb, objective_names = keys(obj_trajs), trajectories = obj_trajs))
        obs_attitude = compute_attitude(spec, data[obs_name], sun_dir)

        # Orientación de la cámara respecto a ECI en cada instante: C_ci = C_cb · C_bi
        C_ci = camera_dcm_history(C_cb, obs_attitude)
        camera_sim!(data, los_data, photometry_data, C_ci, obs_name, env.obj_keys, output_path, env.camera_config)
    end

    export_data_csvs!(data, los_data, env; output_root=output_path)

    println("Execution time: ", time() - t0, " s")
end

"""
    run_spec(spec_path, output_path) -> Nothing

Ejecuta la simulación completa a partir de una spec: genera el JSON completo en
`config/generated/<nombre de la spec>.json` y simula con él.

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
include("scripts/main.jl")
run_spec("config/specs/ejemplo_actitud.json")
```
"""
function run_spec(spec_path::AbstractString, output_path::Union{AbstractString,Nothing} = nothing)
    spec_path = abspath(spec_path)
    name = splitext(basename(spec_path))[1]                # nombre de la spec, sin .json

    # Sin carpeta indicada: una por escenario, para no mezclar resultados de specs distintas
    output_path === nothing && (output_path = joinpath(@__DIR__, "..", "output", name))

    # Un JSON completo tiene los elementos (SMA...) en cada objetivo; una spec tiene "orbit"/"pattern"
    spec = Scenario.read_json(spec_path)
    objs = get(spec, "objectives", Dict{String,Any}())
    any(o -> o isa AbstractDict && haskey(o, "SMA"), values(objs)) && error(
        "'$spec_path' is a full scenario JSON, not a spec. Pass the spec from config/specs/ " *
        "(or call run_simulation(path, output_path) to simulate a full JSON).")

    # 1. Spec -> JSON completo en config/generated/ (mismo nombre que la spec)
    cfg_path = normpath(joinpath(@__DIR__, "..", "config", "generated", name * ".json"))
    generate_config_file(spec_path, cfg_path)

    # 2. Simulación con el JSON generado
    run_simulation(cfg_path, normpath(abspath(output_path)))
    return nothing
end

# Solo se ejecuta al lanzar el script directamente, no al hacer include (VS Code, MonteCarlo)
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    args = parse_arguments()
    run_spec(args.spec_path, args.output_path)
end
