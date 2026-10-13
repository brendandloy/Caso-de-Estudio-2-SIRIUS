push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Project
using Project.Types
using Project.Propagator
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
        attitude_data[obj_name] = propagate_attitude(env.obj_w[obj_name], Time_samples)
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

function run_simulation()
    args = parse_arguments(String[])
    return run_simulation(
        args.config_path,
        args.output_path
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    args = parse_arguments()

    run_simulation(
        args.config_path,
        args.output_path
    )
end