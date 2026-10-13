# =============================================================================
# Kalman input preparation script
# Self-contained script: imports, functions and main
# =============================================================================

using CSV
using DataFrames
using StaticArrays
using LinearAlgebra
using Printf
using DifferentialEquations
using Dates
using GLMakie
set_theme!(theme_light())
using JSON3
using Statistics

# =============================================================================
# DATA STRUCTURES
# =============================================================================

struct StateTrack
    t::Vector{Float64}
    r::Vector{SVector{3,Float64}}
    v::Vector{SVector{3,Float64}}
end

struct LosTrack
    t::Vector{Float64}
    rho::Vector{SVector{3,Float64}}
    Vmag::Vector{Float64}
end

struct LosDatum
    t::Float64
    obs_name::Symbol
    obj_name::Symbol
    rho::SVector{3,Float64}
    Vmag::Float64
end

struct KalmanInputs
    observer_tracks::Dict{Symbol,StateTrack}
    los_tracks::Dict{Tuple{Symbol,Symbol},LosTrack}
    measurements_by_time::Dict{Float64,Vector{LosDatum}}
end

struct ValidationData
    observer_truth::Dict{Symbol,StateTrack}
    rso_truth::Dict{Symbol,StateTrack}
end

struct UKFParameters
    lambda::Float64
    wm::Vector{Float64}
    wc::Vector{Float64}
end

# =============================================================================
# BASIC UTILITIES
# =============================================================================
function arcsec_to_rad(x)
    return x * pi / (180.0 * 3600.0)
end


function parse_float_value(x)
    if x isa Missing
        return NaN
    elseif x isa Number
        return Float64(x)
    else
        s = strip(String(x))
        s = replace(s, "," => ".")
        return parse(Float64, s)
    end
end

function push_state!(track::StateTrack, t::Float64, x::Vector{Float64})
    r = SVector{3,Float64}(x[1], x[2], x[3])
    v = SVector{3,Float64}(x[4], x[5], x[6])

    push!(track.t, t)
    push!(track.r, r)
    push!(track.v, v)
end


function read_semicolon_csv(path::AbstractString)
    isfile(path) || error("File not found: $path")
    return CSV.read(path, DataFrame; delim=';')
end

function require_columns(df::DataFrame, cols::Vector{Symbol}, path::AbstractString)
    missing_cols = setdiff(cols, Symbol.(names(df)))
    isempty(missing_cols) || error("Missing columns in $path: $missing_cols")
    return nothing
end

function vector3_from_row(row, cx::Symbol, cy::Symbol, cz::Symbol)
    return SVector{3,Float64}(
        parse_float_value(row[cx]),
        parse_float_value(row[cy]),
        parse_float_value(row[cz])
    )
end

function normalize_los(v::SVector{3,Float64})
    n = norm(v)
    n > 0.0 || error("LOS vector with zero norm found.")
    return v / n
end

function clean_symbol_name(x)
    return Symbol(strip(String(x)))
end

function format_decimal_2(value::Real)
    return replace(@sprintf("%.2f", Float64(value)), "." => ",")
end

# =============================================================================
# STATE MEASUREMENT PARSER
# =============================================================================

function parse_observer_state_file(path::AbstractString)
    df = read_semicolon_csv(path)

    meas_cols = [
        :t,
        :r_x_meas, :r_y_meas, :r_z_meas,
        :v_x_meas, :v_y_meas, :v_z_meas
    ]

    true_cols = [
        :t,
        :r_x_true, :r_y_true, :r_z_true,
        :v_x_true, :v_y_true, :v_z_true
    ]

    require_columns(df, meas_cols, path)
    require_columns(df, true_cols, path)

    sort!(df, :t)

    t = parse_float_value.(df.t)

    r_meas = SVector{3,Float64}[]
    v_meas = SVector{3,Float64}[]

    r_true = SVector{3,Float64}[]
    v_true = SVector{3,Float64}[]

    for row in eachrow(df)
        push!(r_meas, vector3_from_row(row, :r_x_meas, :r_y_meas, :r_z_meas))
        push!(v_meas, vector3_from_row(row, :v_x_meas, :v_y_meas, :v_z_meas))

        push!(r_true, vector3_from_row(row, :r_x_true, :r_y_true, :r_z_true))
        push!(v_true, vector3_from_row(row, :v_x_true, :v_y_true, :v_z_true))
    end

    measured_track = StateTrack(t, r_meas, v_meas)
    truth_track = StateTrack(t, r_true, v_true)

    return measured_track, truth_track
end

# =============================================================================
# LOS MEASUREMENT PARSER
# =============================================================================

function parse_los_file(
    path::AbstractString,
    rso_name::AbstractString,
    max_vmag::Real
)
    df = read_semicolon_csv(path)

    required_cols = [
        :t, :obj_name,
        :dir_x_meas, :dir_y_meas, :dir_z_meas,
        :Vmag
    ]

    require_columns(df, required_cols, path)

    target_obj_name = "obj_" * strip(rso_name)

    df = df[String.(df.obj_name) .== target_obj_name, :]

    df.Vmag = parse_float_value.(df.Vmag)

    mask_vmag = isfinite.(df.Vmag) .& (df.Vmag .<= max_vmag)
    df = df[mask_vmag, :]

    if nrow(df) == 0
        println(
            "No LOS measurements found for RSO '$target_obj_name' " *
            "with Vmag <= $max_vmag in $path. This observer will be skipped."
        )

        return Dict{Symbol,LosTrack}()
    end

    sort!(df, :t)

    obj_key = Symbol(target_obj_name)

    t = Float64[]
    rho = SVector{3,Float64}[]
    Vmag = Float64[]

    for row in eachrow(df)
        rho_raw = vector3_from_row(
            row,
            :dir_x_meas,
            :dir_y_meas,
            :dir_z_meas
        )

        push!(t, parse_float_value(row[:t]))
        push!(rho, normalize_los(rho_raw))
        push!(Vmag, row[:Vmag])
    end

    return Dict(obj_key => LosTrack(t, rho, Vmag))
end

# =============================================================================
# TRUTH PARSER
# =============================================================================

function parse_rso_truth_file(path::AbstractString, rso_name::AbstractString)
    if !isfile(path)
        @warn "RSO truth file not found. Validation truth will be empty." path
        return Dict{Symbol,StateTrack}()
    end

    df = read_semicolon_csv(path)

    required_cols = [
        :t, :obj_name,
        :r_x, :r_y, :r_z,
        :v_x, :v_y, :v_z
    ]

    require_columns(df, required_cols, path)

    target_obj_name = "obj_" * strip(rso_name)

    mask = String.(df.obj_name) .== target_obj_name
    df = df[mask, :]

    if nrow(df) == 0
        error("No truth states found for RSO '$target_obj_name' in $path")
    end

    sort!(df, :t)

    obj_key = Symbol(target_obj_name)

    t = Float64[]
    r = SVector{3,Float64}[]
    v = SVector{3,Float64}[]

    for row in eachrow(df)
        push!(t, parse_float_value(row[:t]))
        push!(r, vector3_from_row(row, :r_x, :r_y, :r_z))
        push!(v, vector3_from_row(row, :v_x, :v_y, :v_z))
    end

    return Dict(obj_key => StateTrack(t, r, v))
end

# =============================================================================
# DIRECTORY DISCOVERY
# =============================================================================

function find_observer_dirs(output_dir::AbstractString, obs_prefix::AbstractString, observer_names)
    obs_dirs = String[]

    for observer_name in observer_names
        obs_dir = joinpath(output_dir, obs_prefix * strip(observer_name))

        isdir(obs_dir) || error("Observer folder not found: $obs_dir")

        push!(obs_dirs, obs_dir)
    end

    return obs_dirs
end

function build_measurement_time_index(los_tracks::Dict{Tuple{Symbol,Symbol},LosTrack})
    by_time = Dict{Float64,Vector{LosDatum}}()

    for ((obs_name, obj_name), track) in los_tracks
        for k in eachindex(track.t)
            datum = LosDatum(
                track.t[k],
                obs_name,
                obj_name,
                track.rho[k],
                track.Vmag[k]
            )

            push!(get!(by_time, track.t[k], LosDatum[]), datum)
        end
    end

    return by_time
end

# =============================================================================
# INPUT LOADER
# =============================================================================

function load_kalman_inputs(cfg, input_dir::AbstractString)

    observer_tracks = Dict{Symbol,StateTrack}()
    observer_truth = Dict{Symbol,StateTrack}()
    los_tracks = Dict{Tuple{Symbol,Symbol},LosTrack}()

    obs_dirs = find_observer_dirs(input_dir, cfg.obs_dir_prefix, cfg.filter_observers)

    for obs_dir in obs_dirs
        obs_name = Symbol(basename(obs_dir))

        state_path = joinpath(obs_dir, cfg.state_file)
        los_path = joinpath(obs_dir, cfg.los_file)

        measured_state, truth_state = parse_observer_state_file(state_path)
        observer_tracks[obs_name] = measured_state

        if cfg.load_truth
            observer_truth[obs_name] = truth_state
        end

        obj_los_tracks = parse_los_file(los_path, cfg.target_rso, cfg.max_vmag)

        for (obj_name, track) in obj_los_tracks
            los_tracks[(obs_name, obj_name)] = track
        end
    end

    measurements_by_time = build_measurement_time_index(los_tracks)

    rso_truth = cfg.load_truth ? parse_rso_truth_file(joinpath(input_dir, cfg.rso_truth_file), cfg.target_rso) :
                             Dict{Symbol,StateTrack}()

    inputs = KalmanInputs(
        observer_tracks,
        los_tracks,
        measurements_by_time
    )

    truth = ValidationData(
        observer_truth,
        rso_truth
    )

    return inputs, truth
end

# =============================================================================
# KALMAN FILTER
# =============================================================================

##-------------------INITIALIZATION-------------------
function get_observer_position_at_time(inputs::KalmanInputs, obs_name::Symbol, t::Float64)
    track = inputs.observer_tracks[obs_name]

    idx = findfirst(ti -> ti == t, track.t)

    idx === nothing && error("Epoch $t not found for observer $obs_name.")

    return track.r[idx]
end


function triangulate_rso_position(inputs::KalmanInputs, measurements, t::Float64, cfg)
    length(measurements) >= 2 || error("At least 2 measurements are required for triangulation.")

    n_lines = length(measurements)

    G = zeros(3 * n_lines, 3 + n_lines)
    d = zeros(3 * n_lines)

    for (i, m) in enumerate(measurements)
        r_obs = get_observer_position_at_time(inputs, m.obs_name, t)
        u = m.rho

        row = 3 * (i - 1)

        G[row + 1, 1] = 1.0
        G[row + 2, 2] = 1.0
        G[row + 3, 3] = 1.0

        G[row + 1, 3 + i] = -u[1]
        G[row + 2, 3 + i] = -u[2]
        G[row + 3, 3 + i] = -u[3]

        d[row + 1] = r_obs[1]
        d[row + 2] = r_obs[2]
        d[row + 3] = r_obs[3]
    end

    sol = pinv(G) * d

    return SVector{3,Float64}(sol[1], sol[2], sol[3])
end


function build_initial_covariance(cfg)
    sigma_r = Float64.(cfg.initial_position_sigma_km)
    sigma_v = Float64.(cfg.initial_velocity_sigma_km_s)

    sigmas = vcat(sigma_r, sigma_v)

    return Matrix(Diagonal(sigmas .^ 2))
end


# function init_state_guess(inputs::KalmanInputs, t::Float64, cfg)
#     measurements0 = inputs.measurements_by_time[t]

#     r0 = triangulate_rso_position(inputs, measurements0, t, cfg)

#     measurement_epochs = sort(collect(keys(inputs.measurements_by_time)))

#     t1 = nothing
#     r1 = nothing

#     for tj in measurement_epochs
#         if tj <= t
#             continue
#         end

#         measurements_j = inputs.measurements_by_time[tj]

#         if length(measurements_j) >= 2
#             t1 = tj
#             r1 = triangulate_rso_position(inputs, measurements_j, tj, cfg)
#             break
#         end
#     end

#     if t1 === nothing
#         error("Cannot initialize UKF: no later epoch with at least 2 measurements to estimate initial velocity.")
#     end

#     dt = t1 - t

#     if dt <= 0.0
#         error("Cannot initialize UKF: invalid dt for initial velocity.")
#     end

#     v0 = (r1 - r0) / dt

#     println("Initial velocity estimated using dt = ", dt, " s.")

#     x0 = vcat(collect(r0), collect(v0))

#     P0 = build_initial_covariance(cfg)

#     return x0, P0
# end

function init_state_guess(inputs::KalmanInputs, t::Float64, cfg)

    # Initial position
    measurements0 = inputs.measurements_by_time[t]
    r0 = triangulate_rso_position(inputs, measurements0, t, cfg)

    measurement_epochs = sort(collect(keys(inputs.measurements_by_time)))

    valid_times = Float64[]
    valid_positions = SVector{3,Float64}[]

    # Find the next two epochs with at least two measurements
    for tj in measurement_epochs
        tj <= t && continue

        measurements_j = inputs.measurements_by_time[tj]

        if length(measurements_j) >= 2
            push!(valid_times, tj)
            push!(valid_positions,
                triangulate_rso_position(inputs, measurements_j, tj, cfg)
            )

            length(valid_times) == 2 && break
        end
    end

    length(valid_times) == 2 ||
        error("Cannot initialize UKF: two later epochs with at least 2 measurements are required.")

    t1, t2 = valid_times
    r1, r2 = valid_positions

    h1 = t1 - t
    h2 = t2 - t

    # Quadratic interpolation derivative evaluated at t0
    v0 =
        -(h1 + h2) / (h1 * h2) * r0 +
         h2 / (h1 * (h2 - h1)) * r1 -
         h1 / (h2 * (h2 - h1)) * r2

    println(
        "Initial velocity estimated using quadratic approximation: ",
        "dt1 = ", h1, " s, dt2 = ", h2, " s."
    )

    x0 = vcat(collect(r0), collect(v0))
    P0 = build_initial_covariance(cfg)

    return x0, P0
end


function build_process_noise_covariance(cfg, dt::Float64)
    if !cfg.process_noise_enabled
        return zeros(6, 6)
    end

    # Standard deviation of the unmodeled acceleration in each ECI axis [km/s²]
    sigma_a = Float64.(cfg.process_accel_sigma_km_s2)

    Q = zeros(6, 6)

    for i in 1:3
        # Acceleration variance [km²/s⁴]
        q = sigma_a[i]^2

        # Position-position covariance [km²]
        Q[i, i] = q * dt^4 / 4.0

        # Position-velocity cross-covariance [km²/s]
        Q[i, i + 3] = q * dt^3 / 2.0
        Q[i + 3, i] = q * dt^3 / 2.0

        # Velocity-velocity covariance [km²/s²]
        Q[i + 3, i + 3] = q * dt^2
    end

    return Q
end

# # Remote sensig process noise covariance
# function build_process_noise_covariance(cfg, dt::Float64)
#     if !cfg.process_noise_enabled
#         return zeros(6, 6)
#     end

#     a_pert = Float64.(cfg.process_accel_sigma_km_s2)

#     sigma_r = 0.5 .* a_pert .* dt^2
#     sigma_v = a_pert .* dt

#     Q = zeros(6, 6)

#     Q[1, 1] = sigma_r[1]
#     Q[2, 2] = sigma_r[2]
#     Q[3, 3] = sigma_r[3]

#     Q[4, 4] = sigma_v[1]
#     Q[5, 5] = sigma_v[2]
#     Q[6, 6] = sigma_v[3]

#     return Q
# end


##-------------------A PRIORI STEP-------------------

function build_ukf_parameters(cfg, n::Int)
    alpha = cfg.ukf_alpha
    beta = cfg.ukf_beta
    kappa = cfg.ukf_kappa

    lambda = alpha^2 * (n + kappa) - n

    wm = fill(1.0 / (2.0 * (n + lambda)), 2 * n + 1)
    wc = fill(1.0 / (2.0 * (n + lambda)), 2 * n + 1)

    wm[1] = lambda / (n + lambda)
    wc[1] = lambda / (n + lambda) + 1.0 - alpha^2 + beta

    return UKFParameters(lambda, wm, wc)
end


function make_positive_definite(P::Matrix{Float64}, cfg)
    P_sym = Symmetric(0.5 * (P + P'))

    jitter = cfg.covariance_jitter

    for _ in 1:8
        try
            cholesky(P_sym + jitter * I)
            return Matrix(P_sym + jitter * I)
        catch
            jitter *= 10.0
        end
    end

    error("Covariance matrix is not positive definite.")
end


function generate_sigma_points(x::Vector{Float64}, P::Matrix{Float64}, ukf_params::UKFParameters, cfg)
    n = length(x)

    P_safe = make_positive_definite(P, cfg)

    S = cholesky(Symmetric((n + ukf_params.lambda) * P_safe)).L

    sigma_points = Vector{Vector{Float64}}(undef, 2 * n + 1)

    sigma_points[1] = copy(x)

    for i in 1:n
        delta = S[:, i]

        sigma_points[i + 1] = x + delta
        sigma_points[i + 1 + n] = x - delta
    end

    return sigma_points
end


function load_coefs_n4()

    # Coeficientes EGM2008 desnormalizados hasta grado/orden 4
    C = zeros(5, 5)
    S = zeros(5, 5)

    # --- Coeficientes Cnm ---
    C_vals = [
        (2, 0, -1.0826261739e-03),
        (2, 1, -2.6673947524e-10),  (2, 2,  1.5746153257e-06),
        (3, 0,  2.5324105186e-06),  (3, 1,  2.1931496313e-06),
        (3, 2,  3.0904390039e-07),  (3, 3,  1.0058351341e-07),
        (4, 0,  1.6198975999e-06),  (4, 1, -5.0864356044e-07),
        (4, 2,  7.8374545740e-08),  (4, 3,  5.9215017764e-08),
        (4, 4, -3.9832042487e-09)
    ]

    # --- Coeficientes Snm ---
    S_vals = [
        (2, 0,  0.0000000000e+00),
        (2, 1,  1.7872706485e-09),  (2, 2, -9.0387278920e-07),
        (3, 0,  0.0000000000e+00),  (3, 1,  2.6808708940e-07),
        (3, 2, -2.1143062093e-07),  (3, 3,  1.9722158184e-07),
        (4, 0,  0.0000000000e+00),  (4, 1, -4.4926543214e-07),
        (4, 2,  1.4813503725e-07),  (4, 3, -1.2009461262e-08),
        (4, 4,  6.5246728719e-09)
    ]


   # --- Cargar valores en las matrices ---
    for (n, m, val) in C_vals
        C[n+1, m+1] = val
    end
    for (n, m, val) in S_vals
        S[n+1, m+1] = val
    end

    return C, S
end


function accel_earth_perturbations(r_eci::SVector{3,Float64},
                                   t::Float64,
                                   C, S, θ₀; Nmax::Int = 4)

"""
Calcula la aceleración gravitatoria de la Tierra en coordenadas ECI,
usando las ecuaciones (3.27)–(3.33) de Montenbruck & Gill (2000).

# Parámetros
- `r_ecef`: posición en el marco Tierra-fijo [m].
- `Nmax`: grado máximo del modelo armónico (por defecto 4).
- `θ_g`: ángulo sideral de rotación de la Tierra [rad] (para rotar a ECI).

# Devuelve
- `SVector{3,Float64}` aceleración en ECI [m/s²].
"""

    # --- Rotar ECI → ECEF ---

    U = Rz(θ_g(θ₀, t))
    r_ecef = U * r_eci

    x, y, z = r_ecef
    r2 = x^2 + y^2 + z^2
    r  = sqrt(r2)

    # Factores que se repiten en las recursiones
    xR_r2 = x * R_e / r2
    yR_r2 = y * R_e / r2
    zR_r2 = z * R_e / r2
    R2_r2 = (R_e^2) / r2

    # --- Vnm y Wnm hasta n = Nmax+1 ---
    # Dimensión: índice (n,m) -> [n+1, m+1]
    V = zeros(Float64, Nmax+2, Nmax+2)
    W = zeros(Float64, Nmax+2, Nmax+2)

    # Condiciones iniciales (3.31): V00 = R/r, W00 = 0
    V[1,1] = R_e / r
    W[1,1] = 0.0

    # 1) ZONALES Vn0 con la 3.30 (m = 0), Wn0 = 0
    #    Ojo con n=1 porque aparece V_{-1,0}; se pone a 0
    for n in 1:(Nmax+1)
        if n == 1
            # De 3.30 con V_{-1,0}=0 → simplifica mucho
            V[2,1] = (2*1 - 1)/(1 - 0) * zR_r2 * V[1,1]
        else
            num1 = (2n - 1) / (n - 0)
            num2 = (n + 0 - 1) / (n - 0)
            V[n+1,1] = num1 * zR_r2 * V[n,1] - num2 * R2_r2 * V[n-1,1]
        end
        W[n+1,1] = 0.0
    end

    # 2) Para cada m ≥ 1:
    #    - primero diagonal Vmm, Wmm (3.29)
    #    - luego columna n > m con 3.30
    for m in 1:(Nmax+1)
        # Diagonal Vmm, Wmm (3.29)
        V[m+1, m+1] = (2m - 1) * ( xR_r2 * V[m,m] - yR_r2 * W[m,m] )
        W[m+1, m+1] = (2m - 1) * ( xR_r2 * W[m,m] + yR_r2 * V[m,m] )

        # Columna n = m+1 ... Nmax+1   (3.30)
        for n in (m+1):(Nmax+1)
            num1 = (2n - 1) / (n - m)
            num2 = (n + m - 1) / (n - m)

            V[n+1, m+1] = num1 * zR_r2 * V[n,   m+1] - num2 * R2_r2 * V[n-1, m+1]
            W[n+1, m+1] = num1 * zR_r2 * W[n,   m+1] - num2 * R2_r2 * W[n-1, m+1]
        end
    end

    # --- Aceleraciones en ECEF (3.15–3.17 / 3.18–3.19) ---
    gm_R2 = μ / (R_e^2)

    ax = 0.0
    ay = 0.0
    az = 0.0

    for n in 2:Nmax
        for m in 0:n
            Cnm = C[n+1, m+1]
            Snm = S[n+1, m+1]

            if m == 0
                # Caso m = 0
                ax += gm_R2 * (-Cnm * V[n+2, 2])  # V_{n+1,1}
                ay += gm_R2 * (-Cnm * W[n+2, 2])  # W_{n+1,1}
            else
                fac = (n - m + 2) * (n - m + 1)   # (n-m+2)! / (n-m)!

                ax += gm_R2 * 0.5 * (
                    -Cnm * V[n+2, m+2] - Snm * W[n+2, m+2] +         # V_{n+1,m+1}
                    fac * ( Cnm * V[n+2, m]   + Snm * W[n+2, m] )    # V_{n+1,m-1}
                )

                ay += gm_R2 * 0.5 * (
                    -Cnm * W[n+2, m+2] + Snm * V[n+2, m+2] +
                    fac * (-Cnm * W[n+2, m]   + Snm * V[n+2, m] )
                )
            end

            # z siempre igual, m = 0 incluido
            az += gm_R2 * (n - m + 1) * (
                -Cnm * V[n+2, m+1] - Snm * W[n+2, m+1]   # V_{n+1,m}
            )
        end
    end

    a_ecef = SVector(ax, ay, az)

    # --- Volver a ECI ---
    a_eci = transpose(U) * a_ecef
    return a_eci
end


function theta0_from_utc(jd::Float64)
    era = 2π*(0.7790572732640 + 1.00273781191135448*(jd - jd_J2000))
    return mod(era, 2π)
end

@inline function θ_g(θ₀::Float64, t::Float64)
    return θ₀ + ω_earth * t
end

@inline function Rz(θ::Float64)
    c, s = cos(θ), sin(θ)
    @SMatrix [  c   s   0.0;
               -s   c   0.0;
                0.0 0.0 1.0 ]
end

function propagate_state(x::Vector{Float64}, t_current::Float64, dt::Float64, prop_model)
    r0 = SVector{3,Float64}(x[1], x[2], x[3])
    v0 = SVector{3,Float64}(x[4], x[5], x[6])

    u0 = vcat(collect(r0), collect(v0))

    C = prop_model.C
    S = prop_model.S
    θ₀ = prop_model.theta0
    Nmax = prop_model.Nmax

    function dynamics(u, p, τ)
        r = SVector{3,Float64}(u[1], u[2], u[3])
        v = SVector{3,Float64}(u[4], u[5], u[6])

        t_abs = t_current + τ

        a_central = -μ / norm(r)^3 * r
        # a_j = accel_earth_perturbations(r, t_abs, C, S, θ₀; Nmax=Nmax)

        return vcat(collect(v), collect(a_central)) #+ a_j))
    end

    # Definición del problema
    prob = ODEProblem(dynamics, u0, (0.0, dt))
    sol = solve(prob, Vern9();# Integradores: Tsit5 Vern9 Vern7
        reltol = 1e-11,                      # Relative error tolerance
        abstol = 1e-11,                      # Absolute error tolerance
        dt =  min(60.0, dt),                           # Initial step size guess
        dtmin = 1e-3,                        # Minimum allowable step size
        dtmax = min(2700.0, dt),                      # Maximum allowable step size
        maxiters = 1e7,                      # Maximum number of internal solver iterations
        adaptive = true,                     # Enables adaptive step-size control
        dense = false,                       # Enables continuous solution interpolation
        saveat = [dt],                     # Fixed output sampling interval (independent of internal steps)
        save_everystep = false               # Prevent storing all internal solver steps
    )

    u_final = sol.u[end]

    return Vector{Float64}(u_final)
end


function ukf_predict(x, P, t_current, dt, Q, ukf_params, cfg, prop_model)

    sigma_points = generate_sigma_points(x, P, ukf_params, cfg)

    propagated_sigma_points = [
        propagate_state(sp, t_current, dt, prop_model) for sp in sigma_points
    ]

    n = length(x)

    x_pred = zeros(n)

    for i in eachindex(propagated_sigma_points)
        x_pred += ukf_params.wm[i] * propagated_sigma_points[i]
    end

    P_pred = zeros(n, n)

    for i in eachindex(propagated_sigma_points)
        dx = propagated_sigma_points[i] - x_pred
        P_pred += ukf_params.wc[i] * (dx * dx')
    end

    P_pred += Matrix(Q)
    P_pred = 0.5 * (P_pred + P_pred')

    return x_pred, P_pred, propagated_sigma_points
end


##-------------------MEASUREMENTS-------------------

function measurement_model(x::Vector{Float64}, inputs::KalmanInputs, measurements, t::Float64)
    if measurements === nothing || isempty(measurements)
        return Float64[]
    end

    r_rso = SVector{3,Float64}(x[1], x[2], x[3])

    y_pred = Float64[]

    for m in measurements
        r_obs = get_observer_position_at_time(inputs, m.obs_name, t)

        rho_vec = r_rso - r_obs
        rho = rho_vec / norm(rho_vec)

        append!(y_pred, collect(rho))
    end

    return y_pred
end


function build_measurement_vector(measurements)
    if measurements === nothing || isempty(measurements)
        return Float64[]
    end

    y_meas = Float64[]

    for m in measurements
        append!(y_meas, collect(m.rho))
    end

    return y_meas
end


##-------------------COVARIANCE MATRIXES-------------------


function build_single_los_covariance(cfg, r_rso::SVector{3,Float64}, r_obs::SVector{3,Float64})
    rho_vec = r_rso - r_obs
    range_km = norm(rho_vec)
    rho = rho_vec / range_km

    I3 = Matrix{Float64}(I, 3, 3)

    P_perp = I3 - collect(rho) * collect(rho)'

    sigma_los = arcsec_to_rad(Float64(cfg.los_sigma_arcsec))
    R_los = sigma_los^2 * P_perp

    sigma_obs = Float64.(cfg.observer_position_sigma_km)
    P_obs = Matrix(Diagonal(sigma_obs .^ 2))

    H_obs = -P_perp / range_km
    R_obs = H_obs * P_obs * H_obs'

    R_block = R_los + R_obs
    # Add jitter to ensure positive definiteness
    R_block += cfg.covariance_jitter * I3

    return 0.5 * (R_block + R_block')
end


function build_measurement_noise_covariance(cfg, x_pred::Vector{Float64}, inputs::KalmanInputs, measurements, t::Float64)
    if measurements === nothing || isempty(measurements)
        return zeros(0, 0)
    end

    r_rso = SVector{3,Float64}(x_pred[1], x_pred[2], x_pred[3])

    n_meas = length(measurements)

    R = zeros(3 * n_meas, 3 * n_meas)

    for (i, m) in enumerate(measurements)
        r_obs = get_observer_position_at_time(inputs, m.obs_name, t)

        R_block = build_single_los_covariance(cfg, r_rso, r_obs)

        idx = 3 * (i - 1) + 1 : 3 * i

        R[idx, idx] .= R_block
    end

    return R
end


function build_measurement_update_matrices(x_pred::Vector{Float64}, sigma_points_pred, ukf_params::UKFParameters, inputs::KalmanInputs, measurements, t::Float64, R::AbstractMatrix{Float64})
    if measurements === nothing || isempty(measurements)
        return Float64[], zeros(0, 0), zeros(length(x_pred), 0)
    end

    y_sigma = [
        measurement_model(sp, inputs, measurements, t) for sp in sigma_points_pred
    ]

    y_pred = zeros(length(y_sigma[1]))

    for i in eachindex(y_sigma)
        y_pred += ukf_params.wm[i] * y_sigma[i]
    end

    n_x = length(x_pred)
    n_y = length(y_pred)

    Pyy = zeros(n_y, n_y)
    Pxy = zeros(n_x, n_y)

    for i in eachindex(sigma_points_pred)
        dx = sigma_points_pred[i] - x_pred
        dy = y_sigma[i] - y_pred

        Pyy += ukf_params.wc[i] * (dy * dy')
        Pxy += ukf_params.wc[i] * (dx * dy')
    end

    Pyy += Matrix(R)

    # Enforce symmetry for numerical stability
    Pyy = 0.5 * (Pyy + Pyy')

    return y_pred, Pyy, Pxy
end


function ukf_update(x_pred::Vector{Float64}, P_pred::Matrix{Float64}, y_meas::Vector{Float64}, y_pred::Vector{Float64}, Pyy::Matrix{Float64}, Pxy::Matrix{Float64})
    if isempty(y_meas)
        return x_pred, P_pred, Float64[], zeros(0, 0)
    end

    innovation = y_meas - y_pred

    K = Pxy / Pyy

    x_est = x_pred + K * innovation

    P_est = P_pred - K * Pyy * K'
    P_est = 0.5 * (P_est + P_est')

    return x_est, P_est, innovation, Pyy
end


# =============================================================================
# DIAGNOSTICS
# =============================================================================

function print_inputs_summary(inputs::KalmanInputs, truth::ValidationData)
    println("\n================ KALMAN INPUT SUMMARY ================")

    println("\nObserver measured tracks:")
    for (obs_name, track) in sort(collect(inputs.observer_tracks); by=x -> string(x[1]))
        @printf("  %-12s  epochs = %d  t = [%.3f, %.3f] s\n",
                string(obs_name), length(track.t), first(track.t), last(track.t))
    end

    println("\nLOS measurement tracks:")
    for ((obs_name, obj_name), track) in sort(collect(inputs.los_tracks); by=x -> string(x[1]))
        @printf("  %-12s -> %-12s  measurements = %d  t = [%.3f, %.3f] s\n",
                string(obs_name), string(obj_name), length(track.t), first(track.t), last(track.t))
    end

    println("\nMeasurement epochs:")
    epochs = sort(collect(keys(inputs.measurements_by_time)))
    @printf("  total epochs with LOS = %d\n", length(epochs))
    if !isempty(epochs)
        @printf("  first epoch = %.3f s, last epoch = %.3f s\n", first(epochs), last(epochs))
    end

    println("\nTruth data loaded only for validation:")
    @printf("  observer truth tracks = %d\n", length(truth.observer_truth))
    @printf("  RSO truth tracks      = %d\n", length(truth.rso_truth))

    println("======================================================\n")
end

# =============================================================================
# GRAPHIC DISPLAY AND ERROR METRICS
# =============================================================================


function compute_state_component_errors(est::StateTrack, truth::StateTrack)
    valid_est_idx = Int[]
    valid_truth_idx = Int[]

    for i in eachindex(est.t)
        j = findfirst(ti -> isapprox(ti, est.t[i]; atol=1e-9, rtol=0.0), truth.t)

        if j !== nothing
            push!(valid_est_idx, i)
            push!(valid_truth_idx, j)
        end
    end

    isempty(valid_est_idx) && error("No matching epochs found between estimated and truth tracks.")

    n = length(valid_est_idx)

    t = zeros(n)
    err_r = zeros(n, 3)
    err_v = zeros(n, 3)

    for k in 1:n
        i = valid_est_idx[k]
        j = valid_truth_idx[k]

        t[k] = est.t[i]

        err_r[k, :] .= collect(est.r[i] - truth.r[j])
        err_v[k, :] .= collect(est.v[i] - truth.v[j])
    end

    return t, err_r, err_v, valid_est_idx
end


function extract_component_sigmas(P_hist, valid_idx)
    n = length(valid_idx)

    sigma_r = zeros(n, 3)
    sigma_v = zeros(n, 3)

    for k in 1:n
        i = valid_idx[k]

        d = diag(P_hist[i])

        sigma_r[k, :] .= sqrt.(max.(d[1:3], 0.0))
        sigma_v[k, :] .= sqrt.(max.(d[4:6], 0.0))
    end

    return sigma_r, sigma_v
end


function rotate_dop_geometry(dr_list; angle_deg::Float64 = 0.1)

    θ = deg2rad(angle_deg)

    # Small rotation around x to avoid the azimuth/elevation singularity
    Rx = @SMatrix [
        1.0  0.0      0.0;
        0.0  cos(θ)  -sin(θ);
        0.0  sin(θ)   cos(θ)
    ]

    return [Rx * dr for dr in dr_list]
end


function compute_dop(inputs::KalmanInputs,
                     measurements,
                     t::Float64,
                     r_rso,
                     cfg)

    measurements === nothing && return NaN
    length(measurements) >= 2 || return NaN

    dr_list = SVector{3,Float64}[]

    # Build observer-to-RSO relative position vectors
    for m in measurements
        r_obs = get_observer_position_at_time(inputs, m.obs_name, t)
        dr = r_rso - r_obs

        norm(dr) > 0.0 || continue
        push!(dr_list, dr)
    end

    length(dr_list) >= 2 || return NaN

    # Check azimuth/elevation singularity
    singular = any(dr -> hypot(dr[1], dr[2]) / norm(dr) < 1e-10, dr_list)

    if singular
        dr_list = rotate_dop_geometry(dr_list)
    end

    σ_los = arcsec_to_rad(Float64(cfg.los_sigma_arcsec))

    J = zeros(3, 3)

    for dr in dr_list
        dx, dy, dz = dr

        ρ² = dx^2 + dy^2
        ρ  = sqrt(ρ²)
        r² = ρ² + dz^2
        r  = sqrt(r²)

        ρ > 0.0 || return NaN

        # Elevation angle
        α = asin(dz / r)

        # Zhang azimuth/elevation Jacobian
        H = [
            -dy / ρ²          -(dx * dz) / (r² * ρ);
             dx / ρ²          -(dy * dz) / (r² * ρ);
             0.0               ρ / r²
        ]

        # Angular covariance consistent with isotropic LOS noise
        σ_α = σ_los
        σ_β = σ_los / cos(α)

        R = [
            σ_β^2  0.0;
            0.0    σ_α^2
        ]

        J .+= H * inv(R) * H'
    end

    λ = eigvals(Symmetric(J))

    λmax = maximum(λ)

    if λmax <= 0.0 || minimum(λ) <= 1e-12 * λmax
        return NaN
    end

    return sqrt(sum(1.0 ./ λ))
end


function compute_module(err_r::AbstractMatrix, err_v::AbstractMatrix)
    size(err_r, 2) == 3 || error("Position error matrix must have three columns.")
    size(err_v, 2) == 3 || error("Velocity error matrix must have three columns.")

    # module is evaluated as the Euclidean magnitude of the Cartesian error vector.
    position_module_km = vec(sqrt.(sum(abs2, err_r; dims=2)))
    velocity_module_km_s = vec(sqrt.(sum(abs2, err_v; dims=2)))
    velocity_module_m_s = 1000.0 .* velocity_module_km_s

    return position_module_km, velocity_module_m_s
end


function find_convergence_time(
    t::AbstractVector,
    one_drms::AbstractVector,
    threshold::Real;
    max_time_s::Real = 300.0
)

    threshold <= 0.0 &&
        return (index=nothing, time_s=NaN, time_min=NaN, value=NaN)

    t_elapsed = Float64.(t) .- Float64(first(t))
    values = Float64.(one_drms)

    # Only samples inside the convergence window
    window_mask = (t_elapsed .<= max_time_s) .& isfinite.(values)

    # Samples above the threshold
    above_mask = window_mask .& (values .> threshold)

    # No sample above threshold -> converged from the beginning
    if !any(above_mask)
        idx = findfirst(window_mask)

        return (
            index = idx,
            time_s = 0.0,
            time_min = 0.0,
            value = values[idx]
        )
    end

    # Last sample above threshold
    idx_above = findlast(above_mask)
    idx_conv = idx_above + 1

    # The following sample must be below threshold and inside the time window
    if idx_conv <= length(values) &&
       t_elapsed[idx_conv] <= max_time_s &&
       isfinite(values[idx_conv]) &&
       values[idx_conv] < threshold

        return (
            index = idx_conv,
            time_s = t_elapsed[idx_conv],
            time_min = t_elapsed[idx_conv] / 60.0,
            value = values[idx_conv]
        )
    end

    return (
        index = nothing,
        time_s = NaN,
        time_min = NaN,
        value = NaN
    )
end

function zoom_start_index(t::AbstractVector, convergence)
    isempty(t) && error("Time vector cannot be empty.")

    if convergence.index !== nothing
        return convergence.index
    end

    # Without a convergence point, remove the first 10% of visibility time.
    t0 = Float64(first(t))
    tf = Float64(last(t))
    t_zoom = t0 + 0.10 * (tf - t0)

    idx = findfirst(ti -> Float64(ti) >= t_zoom, t)
    return idx === nothing ? 1 : idx
end


function save_component_plot(
    t_min,
    err,
    sigma,
    idx,
    convergence;
    quantity::String,
    units::String,
    scale::Float64,
    title_suffix::String="",
    filename::String,
    output_dir::String=".",
    last_measurement_min=nothing
)
    labels = ["x", "y", "z"]

    fig = Figure(size=(900, 700))

    for j in 1:3
        ax = Axis(
            fig[j, 1],
            xlabel = j == 3 ? "Tiempo [min]" : "",
            ylabel = "δ" * quantity * "_" * labels[j] * " [" * units * "]",
            title = uppercasefirst(quantity) * " error por componente " *
                    labels[j] * title_suffix,
            xtickformat = values -> format_decimal_2.(values),
            ytickformat = values -> format_decimal_2.(values),
            xautolimitmargin = (0.0, 0.0),
            yautolimitmargin = (0.0, 0.0)
        )

        lines!(
            ax,
            t_min[idx],
            scale .* err[idx, j],
            label="error"
        )

        lines!(
            ax,
            t_min[idx],
            3.0 * scale .* sigma[idx, j],
            linestyle=:dash,
            label="+3σ"
        )

        lines!(
            ax,
            t_min[idx],
            -3.0 * scale .* sigma[idx, j],
            linestyle=:dash,
            label="-3σ"
        )

        if convergence.index !== nothing
            vlines!(
                ax,
                [convergence.time_min],
                color=:red,
                linewidth=1.0,
                linestyle=:dash,
                label="t_conv = $(format_decimal_2(convergence.time_min)) min"
            )
        end

        if last_measurement_min !== nothing
            vlines!(
                ax,
                [last_measurement_min],
                linewidth=1.0,
                linestyle=:dot,
                label="última medida"
            )
        end

        axislegend(ax, position=:rb)
    end

    save(joinpath(output_dir, filename), fig)

    return fig
end

function plot_component_errors(
    t,
    err_r,
    err_v,
    P_hist,
    valid_idx,
    position_convergence,
    velocity_convergence,
    last_measurement_idx;
    output_dir="."
)
    mkpath(output_dir)

    t_min = (t .- first(t)) ./ 60.0
    sigma_r, sigma_v = extract_component_sigmas(P_hist, valid_idx)

    last_measurement_min = t_min[last_measurement_idx]

    position_zoom_idx = zoom_start_index(
        t[1:last_measurement_idx],
        position_convergence
    )

    velocity_zoom_idx = zoom_start_index(
        t[1:last_measurement_idx],
        velocity_convergence
    )

    pos_title = position_convergence.index === nothing ?
        " - after initial 10%" : " - post-convergencia"

    vel_title = velocity_convergence.index === nothing ?
        " - after initial 10%" : " - post-convergencia"

    # Measurement window
    fig_r = save_component_plot(
        t_min, err_r, sigma_r,
        1:last_measurement_idx,
        position_convergence;
        quantity="posición",
        units="km",
        scale=1.0,
        filename="position_component_errors.png",
        output_dir=output_dir
    )

    fig_v = save_component_plot(
        t_min, err_v, sigma_v,
        1:last_measurement_idx,
        velocity_convergence;
        quantity="velocidad",
        units="km/s",
        scale=1.0,
        filename="velocity_component_errors.png",
        output_dir=output_dir
    )

    # Convergence -> last measurement
    fig_r_zoom = save_component_plot(
        t_min, err_r, sigma_r,
        position_zoom_idx:last_measurement_idx,
        position_convergence;
        quantity="position",
        units="m",
        scale=1000.0,
        title_suffix=pos_title,
        filename="position_component_errors_zoom.png",
        output_dir=output_dir
    )

    fig_v_zoom = save_component_plot(
        t_min, err_v, sigma_v,
        velocity_zoom_idx:last_measurement_idx,
        velocity_convergence;
        quantity="velocidad",
        units="m/s",
        scale=1000.0,
        title_suffix=vel_title,
        filename="velocity_component_errors_zoom.png",
        output_dir=output_dir
    )

    # Convergence -> end of real state
    fig_r_prop = save_component_plot(
        t_min, err_r, sigma_r,
        position_zoom_idx:length(t),
        position_convergence;
        quantity="posición",
        units="m",
        scale=1000.0,
        title_suffix=" - propagación",
        filename="position_component_errors_propagation.png",
        output_dir=output_dir,
        last_measurement_min=last_measurement_min
    )

    fig_v_prop = save_component_plot(
        t_min, err_v, sigma_v,
        velocity_zoom_idx:length(t),
        velocity_convergence;
        quantity="velocidad",
        units="m/s",
        scale=1000.0,
        title_suffix=" - propagación",
        filename="velocity_component_errors_propagation.png",
        output_dir=output_dir,
        last_measurement_min=last_measurement_min
    )

    return fig_r, fig_v, fig_r_zoom, fig_v_zoom, fig_r_prop, fig_v_prop
end


function save_module_plot(
    t_min,
    values,
    idx,
    convergence,
    threshold;
    ylabel::String,
    title::String,
    filename::String,
    scale::Float64=1.0,
    threshold_units::String="",
    output_dir::String=".",
    last_measurement_min=nothing
)
    fig = Figure(size=(900, 500))

    ax = Axis(
        fig[1, 1],
        xlabel="Tiempo [min]",
        ylabel=ylabel,
        title=title,
        yticks=LinearTicks(10),
        xtickformat = values -> format_decimal_2.(values),
        ytickformat = values -> format_decimal_2.(values),
        xautolimitmargin=(0.0, 0.0),
        yautolimitmargin=(0.0, 0.0)
    )

    lines!(
        ax,
        t_min[idx],
        scale .* values[idx],
        label="Error"
    )

    if threshold > 0.0
        threshold_plot = scale * threshold

        hlines!(
            ax,
            [threshold_plot],
            color=:red,
            linewidth=1.0,
            label="Límite = $(format_decimal_2(threshold_plot)) $(threshold_units)"
        )
    end

    if threshold > 0.0 && convergence.index !== nothing
        vlines!(
            ax,
            [convergence.time_min],
            color=:red,
            linewidth=1.0,
            linestyle=:dash,
            label="t_conv = $(format_decimal_2(convergence.time_min)) min"
        )
    end

    if last_measurement_min !== nothing
        vlines!(
            ax,
            [last_measurement_min],
            linewidth=1.0,
            linestyle=:dot,
            label="última medida"
        )
    end

    axislegend(ax, position=:rt)

    save(joinpath(output_dir, filename), fig)

    return fig
end


function plot_module_errors(
    t,
    position_module_km,
    velocity_module_m_s,
    accuracy_threshold_km,
    velocity_threshold_m_s,
    position_convergence,
    velocity_convergence,
    last_measurement_idx;
    output_dir="."
)
    mkpath(output_dir)

    t_min = (t .- first(t)) ./ 60.0
    last_measurement_min = t_min[last_measurement_idx]

    position_zoom_idx = zoom_start_index(
        t[1:last_measurement_idx],
        position_convergence
    )

    velocity_zoom_idx = zoom_start_index(
        t[1:last_measurement_idx],
        velocity_convergence
    )

    pos_title = position_convergence.index === nothing ?
        "after initial 10%" : "post-convergencia"

    vel_title = velocity_convergence.index === nothing ?
        "after initial 10%" : "post-convergencia"

    # Measurement window
    fig_pos = save_module_plot(
        t_min,
        position_module_km,
        1:last_measurement_idx,
        position_convergence,
        accuracy_threshold_km;
        ylabel="Error [km]",
        title="Módulo del Error de Posición",
        filename="position_module.png",
        threshold_units="km",
        output_dir=output_dir
    )

    fig_vel = save_module_plot(
        t_min,
        velocity_module_m_s,
        1:last_measurement_idx,
        velocity_convergence,
        velocity_threshold_m_s;
        ylabel="Error [m/s]",
        title="Módulo del Error de Velocidad",
        filename="velocity_module.png",
        threshold_units="m/s",
        output_dir=output_dir
    )

    # Convergence -> last measurement
    fig_pos_zoom = save_module_plot(
        t_min,
        position_module_km,
        position_zoom_idx:last_measurement_idx,
        position_convergence,
        accuracy_threshold_km;
        ylabel="Error [m]",
        title="Módulo del Error de Posición - " * pos_title,
        filename="position_module_zoom.png",
        scale=1000.0,
        threshold_units="m",
        output_dir=output_dir
    )

    fig_vel_zoom = save_module_plot(
        t_min,
        velocity_module_m_s,
        velocity_zoom_idx:last_measurement_idx,
        velocity_convergence,
        velocity_threshold_m_s;
        ylabel="Error [m/s]",
        title="Módulo del Error de Velocidad - " * vel_title,
        filename="velocity_module_zoom.png",
        threshold_units="m/s",
        output_dir=output_dir
    )

    # Convergence -> end of real state
    fig_pos_prop = save_module_plot(
        t_min,
        position_module_km,
        position_zoom_idx:length(t),
        position_convergence,
        accuracy_threshold_km;
        ylabel="Error [m]",
        title="Módulo del Error de Posición - " * pos_title,
        filename="position_module_propagation.png",
        scale=1000.0,
        threshold_units="m",
        output_dir=output_dir,
        last_measurement_min=last_measurement_min
    )

    fig_vel_prop = save_module_plot(
        t_min,
        velocity_module_m_s,
        velocity_zoom_idx:length(t),
        velocity_convergence,
        velocity_threshold_m_s;
        ylabel="Error [m/s]",
        title="Módulo del Error de Velocidad - " * vel_title,
        filename="velocity_module_propagation.png",
        threshold_units="m/s",
        output_dir=output_dir,
        last_measurement_min=last_measurement_min
    )

    return (
        fig_pos,
        fig_vel,
        fig_pos_zoom,
        fig_vel_zoom,
        fig_pos_prop,
        fig_vel_prop
    )
end

function plot_dop(t, dop_hist; output_dir=".")

    mkpath(output_dir)

    isempty(t) && return nothing
    length(t) == length(dop_hist) ||
        error("Time and DOP vectors have different lengths.")

    t_min = (t .- first(t)) ./ 60.0

    valid = isfinite.(dop_hist)
    mean_dop = any(valid) ? mean(dop_hist[valid]) : NaN

    fig = Figure(size=(900, 500))

    ax = Axis(
        fig[1, 1],
        xlabel="Tiempo [min]",
        ylabel="DOP [km]",
        title="DOP Geométrico",
        xtickformat = values -> format_decimal_2.(values),
        ytickformat = values -> format_decimal_2.(values),
        xautolimitmargin = (0.0, 0.0),
        yautolimitmargin = (0.0, 0.0)
    )

    lines!(ax, t_min, dop_hist, label="DOP")

    if isfinite(mean_dop)
        hlines!(
            ax,
            [mean_dop],
            linestyle=:dash,
            linewidth=1.0,
            label="Media = $(format_decimal_2(mean_dop)) km"
        )
    end

    axislegend(ax, position=:rt)

    save(joinpath(output_dir, "dop.png"), fig)

    return fig
end


function plot_observer_count(inputs::KalmanInputs, epochs; output_dir=".")

    t0 = first(epochs)
    t_min = (Float64.(epochs) .- t0) ./ 60.0

    n_observers = [
        length(get(inputs.measurements_by_time, t, LosDatum[]))
        for t in epochs
    ]

    fig = Figure(size=(900, 400))

    ax = Axis(
        fig[1, 1],
        xlabel = "Tiempo [min]",
        ylabel = "Número de observadores",
        title = "Observadores disponibles",
        yticks = 0:max(length(inputs.observer_tracks), 1),
        xtickformat = values -> format_decimal_2.(values),
        ytickformat = values -> format_decimal_2.(values)
    )

    stairs!(ax, t_min, n_observers)

    xlims!(ax, t_min[1], t_min[end])
    ylims!(ax, 0, max(length(inputs.observer_tracks)+1, 1))

    save(joinpath(output_dir, "observer_count.png"), fig)

    return fig
end


function save_error_csvs(
    t,
    err_r,
    err_v,
    P_hist,
    valid_idx,
    position_module_km,
    velocity_module_m_s,
    position_convergence,
    velocity_convergence,
    last_measurement_idx,
    dop_hist,
    cfg;
    output_dir="."
)
    mkpath(output_dir)

    t_s = t .- first(t)
    t_min = t_s ./ 60.0

    # -------------------------------------------------------------------------
    # Detailed CSVs
    # -------------------------------------------------------------------------
    if cfg.save_detailed_csv

        sigma_r, sigma_v = extract_component_sigmas(P_hist, valid_idx)

        component_df = DataFrame(
            time_s = t_s,
            time_min = t_min,

            position_error_x_km = err_r[:, 1],
            position_error_y_km = err_r[:, 2],
            position_error_z_km = err_r[:, 3],

            position_sigma_x_km = sigma_r[:, 1],
            position_sigma_y_km = sigma_r[:, 2],
            position_sigma_z_km = sigma_r[:, 3],

            velocity_error_x_km_s = err_v[:, 1],
            velocity_error_y_km_s = err_v[:, 2],
            velocity_error_z_km_s = err_v[:, 3],

            velocity_sigma_x_km_s = sigma_v[:, 1],
            velocity_sigma_y_km_s = sigma_v[:, 2],
            velocity_sigma_z_km_s = sigma_v[:, 3]
        )

        position_df = DataFrame(
            time_s = t_s,
            time_min = t_min,
            position_module_km = position_module_km,
            position_module_m = 1000.0 .* position_module_km,
            threshold_km = fill(Float64(cfg.accuracy_threshold), length(t_s))
        )

        velocity_df = DataFrame(
            time_s = t_s,
            time_min = t_min,
            velocity_module_km_s = velocity_module_m_s ./ 1000.0,
            velocity_module_m_s = velocity_module_m_s,
            threshold_m_s = fill(Float64(cfg.velocity_threshold), length(t_s))
        )

        CSV.write(joinpath(output_dir, "state_component_errors.csv"), component_df; delim=';')
        CSV.write(joinpath(output_dir, "position_module.csv"), position_df; delim=';')
        CSV.write(joinpath(output_dir, "velocity_module.csv"), velocity_df; delim=';')
    end

    # -------------------------------------------------------------------------
    # Summary CSV
    # -------------------------------------------------------------------------
    if cfg.save_summary_csv


        valid_dop = filter(isfinite, dop_hist)
        mean_dop = isempty(valid_dop) ? missing : mean(valid_dop)
        @printf("Mean DOP: %.3f km\n", mean_dop)

        if position_convergence.index === nothing
            position_time_s   = missing
            position_time_min = missing
            position_mean     = missing
            position_std      = missing
            position_rms      = missing
        else
            idx = position_convergence.index + 1
            pos_post = position_module_km[idx:last_measurement_idx]

            position_time_s   = position_convergence.time_s
            position_time_min = position_convergence.time_min
            position_mean     = 1000.0 * mean(pos_post)
            position_std      = 1000.0 * std(pos_post)
            position_rms      = 1000.0 * sqrt(mean(abs2, pos_post))
        end

        if velocity_convergence.index === nothing
            velocity_time_s   = missing
            velocity_time_min = missing
            velocity_mean     = missing
            velocity_std      = missing
            velocity_rms      = missing
        else
            idx = velocity_convergence.index + 1
            vel_post = velocity_module_m_s[idx:last_measurement_idx]

            velocity_time_s   = velocity_convergence.time_s
            velocity_time_min = velocity_convergence.time_min
            velocity_mean     = mean(vel_post)
            velocity_std      = std(vel_post)
            velocity_rms      = sqrt(mean(abs2, vel_post))
        end

        convergence_df = DataFrame(
            metric = ["position", "velocity"],

            threshold = [
                Float64(cfg.accuracy_threshold),
                Float64(cfg.velocity_threshold)
            ],

            threshold_unit = ["km", "m/s"],

            convergence_time_s = Union{Missing,Float64}[
                position_time_s,
                velocity_time_s
            ],

            convergence_time_min = Union{Missing,Float64}[
                position_time_min,
                velocity_time_min
            ],

            post_convergence_mean = Union{Missing,Float64}[
                position_mean,
                velocity_mean
            ],

            post_convergence_std = Union{Missing,Float64}[
                position_std,
                velocity_std
            ],

            post_convergence_rms = Union{Missing,Float64}[
                position_rms,
                velocity_rms
            ],

            error_unit = ["m", "m/s"],

            converged_within_5_min = [
                position_convergence.index !== nothing,
                velocity_convergence.index !== nothing
            ],

            mean_dop_km = Union{Missing,Float64}[mean_dop, mean_dop]
        )

        CSV.write(joinpath(output_dir, "summary.csv"), convergence_df; delim=';')
    end

    return nothing
end

function load_kalman_config(path::AbstractString)

    isfile(path) || error("Config file not found: $path")
    cfg = JSON3.read(read(path, String))

    # --- Paths ---
    obs_dir_prefix = cfg.obs_dir_prefix
    state_file = cfg.state_file
    los_file = cfg.los_file
    rso_truth_file = cfg.rso_truth_file

    # --- Input data ---
    load_truth = cfg.load_truth
    target_rso = cfg.target_rso
    max_vmag = cfg.max_vmag
    filter_observers = String.(cfg.filter_observers)
    init_date = DateTime(cfg.init_date)

    # --- UKF parameters ---
    ukf_alpha = cfg.ukf_alpha
    ukf_beta = cfg.ukf_beta
    ukf_kappa = cfg.ukf_kappa

    # --- Initial state uncertainty ---
    initial_position_sigma_km = Float64.(cfg.initial_position_sigma_km)
    initial_velocity_sigma_km_s = Float64.(cfg.initial_velocity_sigma_km_s)

    # --- Observer navigation uncertainty ---
    observer_position_sigma_km = Float64.(cfg.observer_position_sigma_km)
    observer_velocity_sigma_km_s = Float64.(cfg.observer_velocity_sigma_km_s)

    # --- LOS uncertainty ---
    los_sigma_arcsec = cfg.los_sigma_arcsec

    # --- Process noise ---
    process_noise_enabled = cfg.process_noise_enabled
    process_accel_sigma_km_s2 = Float64.(cfg.process_accel_sigma_km_s2)

    # --- Numerical parameters ---
    covariance_jitter = cfg.covariance_jitter

    # --- Convergence thresholds ---
    accuracy_threshold = cfg.accuracy_threshold_km
    velocity_threshold = cfg.velocity_threshold_ms

    # --- output flags ---
    plot_results      = cfg.plot_results
    save_detailed_csv = cfg.save_detailed_csv
    save_summary_csv  = cfg.save_summary_csv

    println("""
    Loaded Kalman config:
    - Target RSO: $(target_rso)
    - Observers: $(filter_observers)
    - Maximum Vmag: $(max_vmag)
    - Initial epoch: $(init_date)
    - Process noise: $(process_noise_enabled)
    """)

    return (
        obs_dir_prefix = obs_dir_prefix,
        state_file = state_file,
        los_file = los_file,
        rso_truth_file = rso_truth_file,

        load_truth = load_truth,
        target_rso = target_rso,
        max_vmag = max_vmag,
        filter_observers = filter_observers,
        init_date = init_date,

        ukf_alpha = ukf_alpha,
        ukf_beta = ukf_beta,
        ukf_kappa = ukf_kappa,

        initial_position_sigma_km = initial_position_sigma_km,
        initial_velocity_sigma_km_s = initial_velocity_sigma_km_s,

        observer_position_sigma_km = observer_position_sigma_km,
        observer_velocity_sigma_km_s = observer_velocity_sigma_km_s,

        los_sigma_arcsec = los_sigma_arcsec,

        process_noise_enabled = process_noise_enabled,
        process_accel_sigma_km_s2 = process_accel_sigma_km_s2,

        covariance_jitter = covariance_jitter,

        accuracy_threshold = accuracy_threshold,
        velocity_threshold = velocity_threshold,

        plot_results = plot_results,
        save_detailed_csv = save_detailed_csv,  
        save_summary_csv = save_summary_csv
    )
end


function parse_arguments_kalman(args = ARGS)

    default_input  = joinpath(@__DIR__, "..", "output")
    default_config = joinpath(@__DIR__, "..", "config", "kalman_config.json")
    default_output = joinpath(@__DIR__, "..", "output")

    if isempty(args)
        return (
            input_path = default_input,
            config_path = default_config,
            output_path = default_output
        )
    end

    length(args) == 3 || error(
        "Usage: julia Kalman_RSO.jl [<input_path> <config_path> <output_path>]"
    )

    return (
        input_path = abspath(args[1]),
        config_path = abspath(args[2]),
        output_path = abspath(args[3])
    )
end


# =============================================================================
# MAIN
# =============================================================================


const μ = 398600.4415
const R_e = 6378.1363
const ω_earth = 7.2921150e-5
const jd_J2000 = 2451545.0

function run_kalman(
    input_path::AbstractString,
    config_path::AbstractString,
    output_path::AbstractString
)

    input_path = abspath(input_path)
    config_path = abspath(config_path)
    output_path = abspath(output_path)

    cfg = load_kalman_config(config_path)

    isdir(input_path) || error("Input directory not found: $input_path")
    mkpath(output_path)

    println("Input directory: $input_path")
    println("Output directory: $output_path")

    inputs, truth = load_kalman_inputs(cfg, input_path)

    # print_inputs_summary(inputs, truth)

    x_pred = StateTrack(Float64[], SVector{3,Float64}[], SVector{3,Float64}[])
    x_est  = StateTrack(Float64[], SVector{3,Float64}[], SVector{3,Float64}[])

    dop_hist = Float64[]

    P_pred_hist = Matrix{Float64}[]
    P_est_hist  = Matrix{Float64}[]

    y_pred_hist = Vector{Float64}[]
    y_meas_hist = Vector{Float64}[]

    innovation_hist = Vector{Float64}[]
    S_hist = Matrix{Float64}[]

    C_prop, S_prop = load_coefs_n4()
    jd0 = datetime2julian(cfg.init_date)
    theta0_prop = theta0_from_utc(jd0)

    prop_model = (
        C = C_prop,
        S = S_prop,
        theta0 = theta0_prop,
        Nmax = 4
    )

    obs_ref = first(values(inputs.observer_tracks))
    all_epochs = obs_ref.t

    t0 = findfirst(t -> begin
        meas = get(inputs.measurements_by_time, t, nothing)
        meas !== nothing && length(meas) >= 2
    end, all_epochs)

    t0 === nothing && error("No epoch with at least 2 measurements found.")

    x_0, P0 = init_state_guess(inputs, all_epochs[t0], cfg)
    ukf_params = build_ukf_parameters(cfg, length(x_0))

    let x = x_0, P = P0, t_prev = all_epochs[t0]

        for (n, t) in enumerate(all_epochs[t0:end])

            measurements = get(inputs.measurements_by_time, t, nothing)

            if n == 1
                x_ = x_0
                P_ = P0
                sigma_points_ = generate_sigma_points(x_, P_, ukf_params, cfg)
            else
                dt = t - t_prev
                Q = build_process_noise_covariance(cfg, dt) # Could be a constant (dt cte)
                x_, P_, sigma_points_ = ukf_predict(x, P, t_prev, dt, Q, ukf_params, cfg, prop_model)
            end

            r_rso_dop = SVector{3,Float64}(x_[1], x_[2], x_[3])

            dop = compute_dop(inputs,  measurements, t, r_rso_dop, cfg)

            push!(dop_hist, dop)

            y_ = measurement_model(x_, inputs, measurements, t)

            y = build_measurement_vector(measurements)

            R = build_measurement_noise_covariance(cfg, x_, inputs, measurements, t)

            y_, Pyy, Pxy = build_measurement_update_matrices(x_, sigma_points_, ukf_params, inputs, measurements, t, R)

            x, P, innovation, S = ukf_update(x_, P_, y, y_, Pyy, Pxy)

            push_state!(x_pred, t, x_)
            push_state!(x_est, t, x)

            push!(P_pred_hist, Matrix(P_))
            push!(P_est_hist, Matrix(P))

            push!(y_pred_hist, y_)
            push!(y_meas_hist, y)

            push!(innovation_hist, innovation)
            push!(S_hist, Matrix(S))
            
            t_prev = t
        end
    end


    target_key = Symbol("obj_" * cfg.target_rso)

    truth_track = truth.rso_truth[target_key]

    t_cmp, err_r, err_v, valid_idx = compute_state_component_errors(x_est, truth_track)

    position_module_km, velocity_module_m_s = compute_module(err_r, err_v)

    # Last epoch with at least one LOS measurement
    last_measurement_time = maximum(keys(inputs.measurements_by_time))

    last_measurement_idx = findlast(t_cmp .<= last_measurement_time)

    last_measurement_idx === nothing &&
        error("No measurement epoch found inside the estimation time span.")

    println(
        "Last measurement epoch: ",
        (t_cmp[last_measurement_idx] - first(t_cmp)) / 60.0,
        " min"
    )

    position_convergence = find_convergence_time(
        t_cmp[1:last_measurement_idx],
        position_module_km[1:last_measurement_idx],
        cfg.accuracy_threshold;
        max_time_s = 300.0
    )

    velocity_convergence = find_convergence_time(
        t_cmp[1:last_measurement_idx],
        velocity_module_m_s[1:last_measurement_idx],
        cfg.velocity_threshold;
        max_time_s = 300.0
    )

    if cfg.accuracy_threshold <= 0.0
        println("Position convergence threshold disabled.")
    elseif position_convergence.index === nothing
        println("Position convergence not found within the first 5 minutes.")
    else
        @printf(
            "Position convergence: %.3f min, module = %.6f km\n",
            position_convergence.time_min,
            position_convergence.value
        )
    end

    if cfg.velocity_threshold <= 0.0
        println("Velocity convergence threshold disabled.")
    elseif velocity_convergence.index === nothing
        println("Velocity convergence not found within the first 5 minutes.")
    else
        @printf(
            "Velocity convergence: %.3f min, module = %.6f m/s\n",
            velocity_convergence.time_min,
            velocity_convergence.value
        )
    end

    if cfg.plot_results
        println("Plotting results...")
        fig_r, fig_v, fig_r_zoom, fig_v_zoom = plot_component_errors(
            t_cmp,
            err_r,
            err_v,
            P_est_hist,
            valid_idx,
            position_convergence,
            velocity_convergence,
            last_measurement_idx;
            output_dir = output_path
        )

        fig_pos_module, fig_vel_module, fig_pos_module_zoom, fig_vel_module_zoom = plot_module_errors(
            t_cmp,
            position_module_km,
            velocity_module_m_s,
            cfg.accuracy_threshold,
            cfg.velocity_threshold,
            position_convergence,
            velocity_convergence,
            last_measurement_idx;
            output_dir = output_path
        )

        plot_dop( x_pred.t, dop_hist; output_dir=output_path )

        plot_observer_count(inputs, all_epochs[t0:end]; output_dir = output_path)
    end

    save_error_csvs(
        t_cmp,
        err_r,
        err_v,
        P_est_hist,
        valid_idx,
        position_module_km,
        velocity_module_m_s,
        position_convergence,
        velocity_convergence,
        last_measurement_idx,
        dop_hist,
        cfg;
        output_dir=output_path
    )   

    return (
    t = t_cmp,
    position_error_m = 1000.0 .* position_module_km,
    velocity_error_m_s = velocity_module_m_s,
    position_convergence = position_convergence,
    velocity_convergence = velocity_convergence
    )
end

function run_kalman()
    args = parse_arguments_kalman()
    return run_kalman(
        args.input_path,
        args.config_path,
        args.output_path
    )
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    args = parse_arguments_kalman()

    run_kalman(
        args.input_path,
        args.config_path,
        args.output_path
    )
end
