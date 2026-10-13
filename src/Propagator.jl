module Propagator

export propagate, propagate_attitude

using ..Constants
using ..Transformations
using ..Types

using DifferentialEquations
using StaticArrays
using LinearAlgebra
using Dates

#------------------------------------FUNCTIONS----------------------------------------------

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
    xR_r2 = x * R / r2
    yR_r2 = y * R / r2
    zR_r2 = z * R / r2
    R2_r2 = (R^2) / r2

    # --- Vnm y Wnm hasta n = Nmax+1 ---
    # Dimensión: índice (n,m) -> [n+1, m+1]
    V = zeros(Float64, Nmax+2, Nmax+2)
    W = zeros(Float64, Nmax+2, Nmax+2)

    # Condiciones iniciales (3.31): V00 = R/r, W00 = 0
    V[1,1] = R / r
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
    gm_R2 = μ / (R^2)

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


function propagate(obj_dict::Dict{Symbol, State}, start_time, sim_time, Δt_out)

    results = Dict{Symbol, TrajectoryData}()

    jd = datetime2julian(start_time)
    θ₀ = theta0_from_utc(jd)

    # Coeficientes (no normalizados)
    C, S =  load_coefs_n4() #load_EGM2008_coeffs(Nmax)

    for (name, state) in obj_dict
            println("Propagating $name...")

            # Estado inicial
            u0 = vcat(state.r, state.v)

            # Ecuaciones del movimiento
            function dynamics(u, p, t) # El integrador necesita 3 argumentos aunque no se usen
                r = SVector{3}(u[1:3])
                v = SVector{3}(u[4:6])

                a_central = -μ / norm(r)^3 * r
                a_j = accel_earth_perturbations(r, t, C, S, θ₀)

                return vcat(v, a_central + a_j)
            end

            # Definición del problema
            prob = ODEProblem(dynamics, u0, (0.0, sim_time))
            sol = solve(prob, Vern9();# Integradores: Tsit5 Vern9 Vern7
                reltol = 1e-11,                      # Relative error tolerance
                abstol = 1e-11,                      # Absolute error tolerance
                dt = 60.0,                           # Initial step size guess
                dtmin = 1e-3,                        # Minimum allowable step size
                dtmax = 2700.0,                      # Maximum allowable step size
                maxiters = 1e7,                      # Maximum number of internal solver iterations
                adaptive = true,                     # Enables adaptive step-size control
                dense = false,                       # Enables continuous solution interpolation
                saveat = Δt_out,                     # Fixed output sampling interval (independent of internal steps)
                save_everystep = false               # Prevent storing all internal solver steps
            ) 

            r_data = Vector{SVector{3, Float64}}()
            v_data = Vector{SVector{3, Float64}}()

            for u in sol.u
                push!(r_data, SVector{3}(u[1:3]))
                push!(v_data, SVector{3}(u[4:6]))
            end

            results[name] = TrajectoryData(r_data, v_data, sol.t)

        end

    return results

end

# body axes initially coincide with ECI axes
const DEFAULT_DCM_B2I = @SMatrix [
    1.0 0.0 0.0;
    0.0 1.0 0.0;
    0.0 0.0 1.0
]

function propagate_attitude(w::Tuple{Float64,SVector{3,Float64}}, tvec::AbstractVector{<:Real};
                            DCM0_b2i::SMatrix{3,3,Float64,9} = DEFAULT_DCM_B2I)::AttitudeData

    w_mod, w_dir = w

    DCM_b2i = Vector{SMatrix{3,3,Float64,9}}(undef, length(tvec))
    DCM_b2i[1] = DCM0_b2i

    if iszero(w_mod) || iszero(norm(w_dir))
        @warn "Zero angular velocity. The object will be considered non-rotating."

        fill!(DCM_b2i, DCM0_b2i)
        return AttitudeData(DCM_b2i)
    end

    dt = Float64(tvec[2] - tvec[1])
    delta_DCM = axis_angle_rotation(w_dir, w_mod * dt)

    @inbounds for k in 2:length(tvec)
        DCM_b2i[k] = delta_DCM * DCM_b2i[k - 1]
    end

    return AttitudeData(DCM_b2i)
end

end