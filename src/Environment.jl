module Environment

export load_config, parse_arguments, nav_err!, load_kernels, gs_positions, sun_dir_eci, export_data_csvs!
export eclipseEarthECI, losECI, sunGlareConeECI

using ..Types
using ..Transformations
using ..Constants

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


function parse_arguments(args = ARGS)

    default_config = joinpath(@__DIR__, "..", "config", "config.json")
    default_output = joinpath(@__DIR__, "..", "output")

    if isempty(args)
        return (
            config_path = default_config,
            output_path = default_output
        )
    end

    length(args) == 2 || error(
        "Usage: julia main.jl [<config_path> <output_path>]"
    )

    return (
        config_path = abspath(args[1]),
        output_path = abspath(args[2])
    )
end


function load_config(path::AbstractString)

    cfg = JSON3.read(path)

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
    obs_dir = Dict{Symbol, SVector{3,Float64}}()
    obs_pos_sigma = Dict{Symbol, SVector{3,Float64}}()
    obs_vel_sigma = Dict{Symbol, SVector{3,Float64}}()
    obj_w = Dict{Symbol, Tuple{Float64,SVector{3,Float64}}}()
    obs_keys = Symbol[]
    obj_keys = Symbol[]

    for (name, dict) in cfg.objectives
        key = Symbol("obj_$name")
        objects[key] = kepler_to_cartesian(dict)
        obj_w[key] = (deg2rad(Float64(dict.w_deg_s)),SVector{3,Float64}(dict.w_dir))
        push!(obj_keys, key)
    end

    for (name, dict) in cfg.observers
        key = Symbol("obs_$name")
        objects[key] = kepler_to_cartesian(dict)
        v = SVector{3,Float64}(dict.obs_dir)
        obs_dir[key] = v / norm(v)
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
        obs_dir                = obs_dir,
        obs_pos_sigma          = obs_pos_sigma,
        obs_vel_sigma          = obs_vel_sigma,
        obs_keys               = obs_keys,  
        obj_keys               = obj_keys,
        obj_w                  = obj_w,
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
        # println("Final position of $key: ", out[key].r[end][1:3])
    end

    return out
end


function sun_dir_eci(t0::Float64, t_sim::Float64)
    et = t0 + t_sim
    state, _ = spkezr("SUN", et, "J2000", "NONE", "EARTH")
    v = state[1:3] / norm(state[1:3])
    return SVector{3,Float64}(v)
end


function eclipseEarthECI(r_sc::SVector{3,Float64}, u_sun::SVector{3,Float64}; return_details::Bool=false)
    # Numerical tolerance
    Tol = 1e-9

    # Shadow axis points away from the Sun
    u_axis = -u_sun

    # Spacecraft coordinates in the shadow frame
    x = dot(r_sc, u_axis)
    rho = norm(r_sc - x * u_axis)

    # Spacecraft is on the sunward side of Earth
    if x <= 0.0
        if return_details
            details = (; u_axis=u_axis, x_axis=x, rho_axis=rho, Re_used=R,
                Rs_used=R_sun, Dsun=D_sun, r_umbra=0.0, r_penum=0.0,
                BlockReason=:sunward_hemisphere)
            return false, "lit", details
        end
        return false, "lit", nothing
    end

    # Conical shadow model
    L_umbra = D_sun * R / (R_sun - R)
    L_penum = D_sun * R / (R_sun + R)

    r_umbra = max(R * (1.0 - x / L_umbra), 0.0)
    r_penum = R * (1.0 + x / L_penum)

    if (x <= L_umbra + Tol) && (r_umbra > 0.0) && (rho <= r_umbra + Tol)
        if return_details
            details = (; u_axis=u_axis, x_axis=x, rho_axis=rho, Re_used=R,
                Rs_used=R_sun, Dsun=D_sun, L_umbra=L_umbra, L_penum=L_penum,
                r_umbra=r_umbra, r_penum=r_penum, BlockReason=:inside_umbra_cone)
            return true, "umbra", details
        end
        return true, "umbra", nothing
    elseif rho <= r_penum + Tol
        if return_details
            details = (; u_axis=u_axis, x_axis=x, rho_axis=rho, Re_used=R,
                Rs_used=R_sun, Dsun=D_sun, L_umbra=L_umbra, L_penum=L_penum,
                r_umbra=r_umbra, r_penum=r_penum, BlockReason=:inside_penumbra_cone)
            return true, "penumbra", details
        end
        return true, "penumbra", nothing
    else
        if return_details
            details = (; u_axis=u_axis, x_axis=x, rho_axis=rho, Re_used=R,
                Rs_used=R_sun, Dsun=D_sun, L_umbra=L_umbra, L_penum=L_penum,
                r_umbra=r_umbra, r_penum=r_penum, BlockReason=:outside_shadow_cones)
            return false, "lit", details
        end
        return false, "lit", nothing
    end
end


function losECI(r1::SVector{3,Float64}, r2::SVector{3,Float64}; return_details::Bool=false)
    # Numerical parameters
    tol = 1e-9
    max_project_depth = 100.0

    r1v = r1
    r2v = r2

    r1n = norm(r1v)
    r2n = norm(r2v)

    projected_r1 = false
    projected_r2 = false
    projdepth_r1 = 0.0
    projdepth_r2 = 0.0

    # Soft projection of endpoint 1 onto the Earth sphere
    depth1 = R - r1n
    if depth1 > tol
        if depth1 <= max_project_depth
            r1v = (R / r1n) * r1v
            r1n = R
            projected_r1 = true
            projdepth_r1 = depth1
        else
            if return_details
                details = (; RadiusNominal=R, RadiusUsed=R, r1_norm=r1n, r2_norm=r2n,
                    Projected_r1=projected_r1, Projected_r2=projected_r2,
                    ProjDepth_r1=projdepth_r1, ProjDepth_r2=projdepth_r2,
                    BlockReason=:endpoint_far_below_surface)
                return false, details
            end
            return false, nothing
        end
    end

    # Soft projection of endpoint 2 onto the Earth sphere
    depth2 = R - r2n
    if depth2 > tol
        if depth2 <= max_project_depth
            r2v = (R / r2n) * r2v
            r2n = R
            projected_r2 = true
            projdepth_r2 = depth2
        else
            if return_details
                details = (; RadiusNominal=R, RadiusUsed=R, r1_norm=r1n, r2_norm=r2n,
                    Projected_r1=projected_r1, Projected_r2=projected_r2,
                    ProjDepth_r1=projdepth_r1, ProjDepth_r2=projdepth_r2,
                    BlockReason=:endpoint_far_below_surface)
                return false, details
            end
            return false, nothing
        end
    end

    # Line-sphere intersection coefficients
    d = r2v - r1v
    a = dot(d, d)

    if a <= tol
        if return_details
            details = (; RadiusNominal=R, RadiusUsed=R, r1_norm=r1n, r2_norm=r2n,
                Projected_r1=projected_r1, Projected_r2=projected_r2,
                ProjDepth_r1=projdepth_r1, ProjDepth_r2=projdepth_r2,
                BlockReason=:coincident_points)
            return false, details
        end
        return false, nothing
    end

    b = 2.0 * dot(r1v, d)
    c = dot(r1v, r1v) - R^2
    disc = b * b - 4.0 * a * c

    # No intersection with the infinite line
    if disc < -abs(b * b + 4.0 * a * abs(c)) * tol
        if return_details
            details = (; RadiusNominal=R, RadiusUsed=R, r1_norm=r1n, r2_norm=r2n,
                Projected_r1=projected_r1, Projected_r2=projected_r2,
                ProjDepth_r1=projdepth_r1, ProjDepth_r2=projdepth_r2,
                BlockReason=:no_line_sphere_intersection)
            return true, details
        end
        return true, nothing
    end

    # Tangency case
    if abs(disc) <= max(1.0, b * b + 4.0 * a * abs(c)) * tol
        s = -b / (2.0 * a)

        if s > tol && s < 1.0 - tol
            if return_details
                details = (; RadiusNominal=R, RadiusUsed=R, r1_norm=r1n, r2_norm=r2n,
                    Projected_r1=projected_r1, Projected_r2=projected_r2,
                    ProjDepth_r1=projdepth_r1, ProjDepth_r2=projdepth_r2,
                    s_intervals=s, BlockReason=:grazing_tangency_inside_segment)
                return false, details
            end
            return false, nothing
        else
            if return_details
                details = (; RadiusNominal=R, RadiusUsed=R, r1_norm=r1n, r2_norm=r2n,
                    Projected_r1=projected_r1, Projected_r2=projected_r2,
                    ProjDepth_r1=projdepth_r1, ProjDepth_r2=projdepth_r2,
                    s_intervals=s, BlockReason=:tangency_at_endpoint_or_outside_segment)
                return true, details
            end
            return true, nothing
        end
    end

    # Two intersections
    sqrtDisc = sqrt(disc)
    s1 = (-b - sqrtDisc) / (2.0 * a)
    s2 = (-b + sqrtDisc) / (2.0 * a)

    smin = min(s1, s2)
    smax = max(s1, s2)

    if (smax > tol) && (smin < 1.0 - tol)
        if return_details
            details = (; RadiusNominal=R, RadiusUsed=R, r1_norm=r1n, r2_norm=r2n,
                Projected_r1=projected_r1, Projected_r2=projected_r2,
                ProjDepth_r1=projdepth_r1, ProjDepth_r2=projdepth_r2,
                s_intervals=(smin, smax), BlockReason=:earth_occlusion)
            return false, details
        end
        return false, nothing
    else
        if return_details
            details = (; RadiusNominal=R, RadiusUsed=R, r1_norm=r1n, r2_norm=r2n,
                Projected_r1=projected_r1, Projected_r2=projected_r2,
                ProjDepth_r1=projdepth_r1, ProjDepth_r2=projdepth_r2,
                s_intervals=(smin, smax),
                BlockReason=:intersections_only_at_endpoints_or_outside_segment)
            return true, details
        end
        return true, nothing
    end
end


function sunGlareConeECI(u_sun::SVector{3,Float64}, r_obs::SVector{3,Float64},
    r_obj::SVector{3,Float64}, angMaxDeg::Real; return_details::Bool=false)

    TolDeg = 1e-6

    # Observer-to-object direction
    v_obj = r_obj - r_obs
    u_obj = v_obj / norm(v_obj)

    # Angular separation between Sun and object directions
    c = clamp(dot(u_sun, u_obj), -1.0, 1.0)
    ang_deg = acosd(c)
    blinded = ang_deg <= (angMaxDeg + TolDeg)

    if return_details
        details = (; u_obs_to_sun=u_sun, u_obs_to_obj=u_obj, dot_uSun_uObj=c,
            ang_deg=ang_deg, angMaxDeg=angMaxDeg, TolDeg=TolDeg)
        return blinded, ang_deg, details
    end

    return blinded, ang_deg, nothing
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