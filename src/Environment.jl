module Environment

export load_config, parse_arguments, nav_err!, load_kernels, gs_positions, sun_dir_eci, export_data_csvs!

using ..Types
using ..Transformations
using ..Constants
using ..Attitude: AttitudeSpec, parse_member_attitude    # lectura del bloque "attitude" de cada miembro
using ..Attitude: camera_mount_block, parse_camera_mount # montaje de la cámara de cada observador (C_cb)
using ..Geometry: ShapeSpec, parse_member_shape          # lectura del bloque "shape" de cada objetivo

using StaticArrays
using SPICE
using Dates
using LinearAlgebra
using JSON3
using CSV
using DataFrames

#------------------------------------FUNCTIONS----------------------------------------------

const _KERNELS_LOADED = Ref(false)

function load_kernels(kernel_dir::AbstractString)
    if _KERNELS_LOADED[]
        return
    end

    furnsh(joinpath(kernel_dir, "de440s.bsp"))
    furnsh(joinpath(kernel_dir, "naif0012.tls.txt"))

    _KERNELS_LOADED[] = true
    return nothing
end


"""
    parse_arguments(args=ARGS) -> NamedTuple

Lee los argumentos de `scripts/main.jl`: la spec es obligatoria y la carpeta de
resultados, opcional.

# Argumentos
- `args`: `[spec.json]` o `[spec.json, output_path]`.

# Devuelve
- `(spec_path, output_path)`: rutas absolutas; sin `output_path`, `nothing` (main.jl usa entonces
  `output/<nombre de la spec>/`).

# Excepciones
- `ErrorException`: si no hay argumentos o hay más de dos (el mensaje indica el uso).

# Ejemplo
```julia
parse_arguments(["config/specs/ejemplo_leo.json"])     # output_path = nothing -> output/ejemplo_leo/
```
"""
function parse_arguments(args = ARGS)
    # La spec es obligatoria: así siempre queda claro qué escenario se simula
    1 <= length(args) <= 2 || error(
        "Usage: julia --project=. scripts/main.jl <spec.json> [<output_path>]"
    )

    return (
        spec_path   = abspath(args[1]),
        output_path = length(args) == 2 ? abspath(args[2]) : nothing
    )
end


function load_config(path::AbstractString)

    # Se lee el texto en memoria: JSON3.read(path) mapea el archivo (mmap) y, mientras la sesión
    # de Julia siga abierta, Windows no deja sobrescribirlo (p. ej. al regenerar la spec)
    cfg = JSON3.read(read(path, String))

    # --- Tiempo ---
    Duration_hours         = cfg.Duration_hours
    Data_time_step_seconds = cfg.Data_time_step_seconds
    init_epoch             = DateTime(cfg.init_epoch)
    camera_config          = cfg.camera_config
    photometry_config      = cfg.photometry_config


    

    # --- Partes del código

    navigation_error_flag = cfg.navigation_error_flag
    graphic_display_flag = cfg.graphic_display_flag 
    ignore_fov_for_measurements_flag = cfg.ignore_fov_for_measurements_flag
    image_generation_flag = cfg.image_generation_flag

    camera_config = merge(
        NamedTuple(camera_config),
        (
            star_csv_flag = Bool(cfg.star_csv_flag),
            navigation_error_flag = Bool(navigation_error_flag),
            ignore_fov_for_measurements_flag = Bool(ignore_fov_for_measurements_flag),
            image_generation_flag = Bool(image_generation_flag)
        )
    )

    # --- Objetos orbitales ---
    objects = Dict{Symbol, State}()
    obs_attitude = Dict{Symbol, AttitudeSpec}()     # actitud del satélite de cada observador
    obs_mount = Dict{Symbol, SMatrix{3,3,Float64,9}}()   # montaje de su cámara, C_cb (cuerpo -> cámara)
    obs_pos_sigma = Dict{Symbol, SVector{3,Float64}}()
    obs_vel_sigma = Dict{Symbol, SVector{3,Float64}}()
    obj_attitude = Dict{Symbol, AttitudeSpec}()     # actitud de cada objetivo (modo + parámetros)
    obj_shape = Dict{Symbol, ShapeSpec}()           # forma de cada objetivo (tipo + medidas + semilla)
    obs_keys = Symbol[]
    obj_keys = Symbol[]

    for (name, dict) in cfg.objectives
        key = Symbol("obj_$name")
        objects[key] = kepler_to_cartesian(dict)
        # Bloque "attitude" obligatorio
        obj_attitude[key] = parse_member_attitude(dict, String(name))
        # Bloque "shape" obligatorio; si aún hay algo por sortear, se sortea aquí
        obj_shape[key] = parse_member_shape(dict, String(name))
        push!(obj_keys, key)
    end

    for (name, dict) in cfg.observers
        key = Symbol("obs_$name")
        objects[key] = kepler_to_cartesian(dict)
        # Bloque "attitude" obligatorio; "camera_mount" opcional (por defecto C_cb = I, con aviso)
        obs_attitude[key] = parse_member_attitude(dict, String(name); role = "Observer")
        obs_mount[key] = parse_camera_mount(camera_mount_block(dict, String(name)), String(name))
        obs_pos_sigma[key] = SVector{3,Float64}(dict.position_sigma_km)
        obs_vel_sigma[key] = SVector{3,Float64}(dict.velocity_sigma_km_s)
        push!(obs_keys, key)
    end

    # --- Estaciones LLH ---
    stations_llh = Dict{Symbol, GroundStationLLH}()

    for (name, s) in cfg.stations_llh
        stations_llh[Symbol(name)] =
            GroundStationLLH(s.lat, s.lon, s.h)
    end

    println("""
    Loaded config:
    - Duration (h): $(Duration_hours)h
    - Step (s): $(Data_time_step_seconds)s
    - Epoch: $(init_epoch)
    - Navigation error: $(navigation_error_flag)
    - Graphics: $(graphic_display_flag)
    - Objects: $(keys(objects))
    - Stations: $(keys(stations_llh))
    """)

    return (
        Duration_hours         = Duration_hours,
        Data_time_step_seconds = Data_time_step_seconds,
        init_epoch             = init_epoch,
        navigation_error_flag  = navigation_error_flag,
        graphic_display_flag   = graphic_display_flag,
        camera_config          = camera_config,
        photometry_config      = photometry_config,
        objects                = objects,
        obs_attitude           = obs_attitude,
        obs_mount              = obs_mount,
        obs_pos_sigma          = obs_pos_sigma,
        obs_vel_sigma          = obs_vel_sigma,
        obs_keys               = obs_keys,  
        obj_keys               = obj_keys,
        obj_attitude           = obj_attitude,
        obj_shape              = obj_shape,
        stations_llh           = stations_llh
    )
end


function nav_err!(objects::Dict{Symbol,State}, obs_keys::Vector{Symbol},
                  obs_pos_sigma::Dict{Symbol,SVector{3,Float64}}, obs_vel_sigma::Dict{Symbol,SVector{3,Float64}})

    for obs_key in obs_keys
        st = objects[obs_key]

        σr = obs_pos_sigma[obs_key]
        σv = obs_vel_sigma[obs_key]

        δr = σr .* SVector{3,Float64}(randn(3))
        δv = σv .* SVector{3,Float64}(randn(3))

        nav_name = Symbol("nav_", String(obs_key))

        objects[nav_name] = State(
            st.r + δr,
            st.v + δv
        )
    end

    return nothing
end


function gs_positions(stations::Dict{Symbol,GroundStationLLH}, start_time, tvec::AbstractVector{<:Real})
    out = Dict{Symbol, TrajectoryData}()
    t_out = collect(Float64.(tvec))

    # Compute θ₀ for the start time
    jd = datetime2julian(start_time)
    θ₀ = theta0_from_utc(jd)

    for (key, st) in stations
        # 1) Convert LLH → ECEF (one time per station)
        lat = deg2rad(st.lat_deg)
        lon = deg2rad(st.lon_deg)
        hkm = st.h_m / 1000.0
        r_ecef = llh_to_ecef(lat, lon, hkm)

        positions = rotate_z.(Ref(r_ecef), θ_g.(θ₀, t_out))

        # Velocities = zero (unless you want real ECI velocity later)
        vzero = SVector{3,Float64}(0.0, 0.0, 0.0)
        velocities = fill(vzero, length(t_out))

        # 3×N matrix (each column is one timestamp)
        out[key] = TrajectoryData(positions, velocities, t_out)
    end

    return out
end


function sun_dir_eci(t0::Float64, t_sim::Float64)
    et = t0 + t_sim
    state, _ = spkezr("SUN", et, "J2000", "NONE", "EARTH")
    v = state[1:3] / norm(state[1:3])
    return SVector{3,Float64}(v)
end


function build_rso_state_df(data, obj_keys)

    df = DataFrame(
        t = Float64[],
        obj_name = String[],
        r_x = Float64[], r_y = Float64[], r_z = Float64[],
        v_x = Float64[], v_y = Float64[], v_z = Float64[]
    )

    for obj_name in obj_keys
        traj = data[obj_name]

        for k in eachindex(traj.t)
            r = traj.r[k]
            v = traj.v[k]

            push!(df, (
                traj.t[k],
                String(obj_name),
                r[1], r[2], r[3],
                v[1], v[2], v[3]
            ))
        end
    end

    return df
end


function build_state_df(data, obs_name; nav_flag::Bool)

    true_traj = data[obs_name]
    nav_name = Symbol("nav_", String(obs_name))

    df = DataFrame(
        t = Float64[],

        r_x_true = Float64[], r_y_true = Float64[], r_z_true = Float64[],
        v_x_true = Float64[], v_y_true = Float64[], v_z_true = Float64[]
    )

    if nav_flag
        df.r_x_meas = Float64[]; df.r_y_meas = Float64[]; df.r_z_meas = Float64[]
        df.v_x_meas = Float64[]; df.v_y_meas = Float64[]; df.v_z_meas = Float64[]
    end

    for k in eachindex(true_traj.t)
        t = true_traj.t[k]
        r = true_traj.r[k]
        v = true_traj.v[k]

        row = (
            t = t,
            r_x_true = r[1], r_y_true = r[2], r_z_true = r[3],
            v_x_true = v[1], v_y_true = v[2], v_z_true = v[3]
        )

        if nav_flag
            rm = data[nav_name].r[k]
            vm = data[nav_name].v[k]

            row = merge(row, (
                r_x_meas = rm[1], r_y_meas = rm[2], r_z_meas = rm[3],
                v_x_meas = vm[1], v_y_meas = vm[2], v_z_meas = vm[3]
            ))
        end

        push!(df, row)
    end

    return df
end


function build_los_df(los_meas, obs_name; nav_flag::Bool)

    df = DataFrame(
        t = Float64[],
        obj_name = String[],

        dir_x_true = Float64[],
        dir_y_true = Float64[],
        dir_z_true = Float64[],

        Vmag = Float64[]
    )

    if nav_flag
        df.dir_x_meas = Float64[]
        df.dir_y_meas = Float64[]
        df.dir_z_meas = Float64[]
    end

    if !haskey(los_meas, obs_name)
        return df
    end

    for (obj_name, measurements) in los_meas[obs_name]
        for m in measurements
            row = (
                t = m.t,
                obj_name = String(obj_name),

                dir_x_true = m.dir_true[1],
                dir_y_true = m.dir_true[2],
                dir_z_true = m.dir_true[3],

                Vmag = m.Vmag
            )

            if nav_flag
                row = merge(row, (
                    dir_x_meas = m.dir_meas[1],
                    dir_y_meas = m.dir_meas[2],
                    dir_z_meas = m.dir_meas[3]
                ))
            end

            push!(df, row)
        end
    end

    return df
end


function export_data_csvs!(data, los_meas, env; output_root=joinpath(@__DIR__, "..", "output"))

    nav_flag = env.navigation_error_flag

    # Common RSO truth file
    rso_df = build_rso_state_df(data, env.obj_keys)
    CSV.write(joinpath(output_root, "rso_states.csv"), rso_df; delim=';')

    for obs_name in env.obs_keys
        outdir = joinpath(output_root, String(obs_name))
        mkpath(outdir)

        state_df = build_state_df(data, obs_name; nav_flag=nav_flag)
        CSV.write(joinpath(outdir, "state_measurements.csv"), state_df; delim=';')

        los_df = build_los_df(los_meas, obs_name; nav_flag=nav_flag)
        CSV.write(joinpath(outdir, "los_measurements.csv"), los_df; delim=';')
    end

    return nothing
end


end # module
