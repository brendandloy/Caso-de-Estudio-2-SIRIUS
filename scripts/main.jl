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
#   2. Simula con ese JSON: propagación -> actitud -> fotometría -> cámara -> CSV.
# Sin <output_path>, los resultados se escriben en output/.
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

using Dates
using SPICE
using StaticArrays

load_kernels(joinpath(@__DIR__, "../data"))

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


    # Example photometric model per tracked object.
    # Replace the category/shape assignment here with your own object definitions if needed.
    debris_models = Dict{Symbol, Dict{Symbol,Any}}()
    for obj_name in env.obj_keys
        debris_models[obj_name] = generate_debris(1, 1; deploy_panels = true, shape_type = 4)
    end

    attitude_data = Dict{Symbol,AttitudeData}()
    photometry_data = Dict{Tuple{Symbol,Symbol}, NamedTuple}()
    los_data = Dict{Symbol, Dict{Symbol, Vector{LOSMeasurement}}}()
        
    for obj_name in env.obj_keys
        # Actitud según el modo del bloque "attitude" (spin, random, fixed, pointing...)
        attitude_data[obj_name] = compute_attitude(env.obj_attitude[obj_name], data[obj_name], sun_dir)
        for obs_name in env.obs_keys
            photometry_data[(obs_name, obj_name)] = compute_apparent_magnitude(
                data[obs_name], data[obj_name], sun_dir,  debris_models[obj_name], env.photometry_config, attitude_data[obj_name].DCM_b2i
            )
        end
    end

    for obs_name in env.obs_keys
        camera_sim!(data, los_data, photometry_data, env.obs_dir[obs_name], obs_name, env.obj_keys, output_path, env.camera_config)
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
- `output_path::AbstractString`: carpeta de resultados (por defecto `output/`).

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
function run_spec(spec_path::AbstractString, output_path::AbstractString = joinpath(@__DIR__, "..", "output"))
    spec_path = abspath(spec_path)

    # Un JSON completo tiene los elementos (SMA...) en cada objetivo; una spec tiene "orbit"/"pattern"
    spec = Scenario.read_json(spec_path)
    objs = get(spec, "objectives", Dict{String,Any}())
    any(o -> o isa AbstractDict && haskey(o, "SMA"), values(objs)) && error(
        "'$spec_path' is a full scenario JSON, not a spec. Pass the spec from config/specs/ " *
        "(or call run_simulation(path, output_path) to simulate a full JSON).")

    # 1. Spec -> JSON completo en config/generated/ (mismo nombre que la spec)
    cfg_path = normpath(joinpath(@__DIR__, "..", "config", "generated", splitext(basename(spec_path))[1] * ".json"))
    generate_config_file(spec_path, cfg_path)

    # 2. Simulación con el JSON generado
    run_simulation(cfg_path, abspath(output_path))
    return nothing
end

# Solo se ejecuta al lanzar el script directamente, no al hacer include (VS Code, MonteCarlo)
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    args = parse_arguments()
    run_spec(args.spec_path, args.output_path)
end
