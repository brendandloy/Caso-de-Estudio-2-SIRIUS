# =============================================================================
# test_attitude.jl
#
# Pruebas del módulo Attitude: convenciones del Tema 6: Attitude Kinematics, 
# compatibilidad con el antiguo propagate_attitude, comportamiento de cada modo,
# montaje de la cámara (C_cb) y orientación de la imagen.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_attitude.jl
#
# Es ligero: solo carga Types, Constants, Transformations, Parameters y Attitude.
# La trayectoria de prueba es una órbita circular kepleriana analítica.
#
# Funciones: circular_traj, axis_angle_rotation, old_propagate_attitude
# =============================================================================

using Test
using LinearAlgebra
using Random
using StaticArrays

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Types)           || include(joinpath(ROOT, "src", "Types.jl"))
isdefined(Main, :Constants)       || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :Transformations) || include(joinpath(ROOT, "src", "Transformations.jl"))
isdefined(Main, :Parameters)      || include(joinpath(ROOT, "src", "Parameters.jl"))
isdefined(Main, :Attitude)        || include(joinpath(ROOT, "src", "Attitude", "Attitude.jl"))

using .Types, .Constants, .Transformations, .Parameters, .Attitude
const A = Attitude                     # acceso corto a las funciones internas

"""
    circular_traj(a_km, inc_deg, t) -> TrajectoryData

Órbita circular kepleriana analítica (RAAN = 0, nodo ascendente en +x ECI).

# Argumentos
- `a_km::Real`: radio de la órbita (km).
- `inc_deg::Real`: inclinación (°).
- `t::AbstractVector`: instantes (s).

# Devuelve
- `TrajectoryData` con r (km) y v (km/s) en ECI.

# Excepciones
- Ninguna.
"""
function circular_traj(a_km::Real, inc_deg::Real, t::AbstractVector)
    n = sqrt(μ / a_km^3)                          # movimiento medio (rad/s)
    ci, si = cosd(inc_deg), sind(inc_deg)
    r = [SVector(a_km * cos(n * tk), a_km * sin(n * tk) * ci, a_km * sin(n * tk) * si) for tk in t]
    v = [SVector(-a_km * n * sin(n * tk), a_km * n * cos(n * tk) * ci, a_km * n * cos(n * tk) * si) for tk in t]
    return TrajectoryData(r, v, collect(Float64.(t)))
end

"""
    axis_angle_rotation(axis, angle) -> SMatrix{3,3,Float64,9}

Matriz de rotación activa de ángulo `angle` alrededor de `axis` (fórmula de
Rodrigues: R = cos θ I + (1 - cos θ) ê êᵀ + sin θ [ê×]).

Por qué existe: es una implementación independiente de Quaternions.jl que sirve
de solución analítica de referencia para el modo `spin` y para
`old_propagate_attitude`. Antes estaba en Transformations.jl, pero solo la usaban
estas pruebas.

# Argumentos
- `axis::SVector{3,Float64}`: eje de giro (unitario).
- `angle::Float64`: ángulo (rad).

# Devuelve
- `SMatrix{3,3,Float64,9}` ortogonal con det = +1.

# Excepciones
- Ninguna (no comprueba que `axis` sea unitario).
"""
function axis_angle_rotation(axis::SVector{3,Float64}, angle::Float64)::SMatrix{3,3,Float64,9}
    ux, uy, uz = axis                                 # componentes del eje unitario
    c, s = cos(angle), sin(angle)
    d = 1.0 - c                                       # término (1 - cos θ) de ê êᵀ
    return @SMatrix [
        c + ux*ux*d       ux*uy*d - uz*s   ux*uz*d + uy*s;
        uy*ux*d + uz*s    c + uy*uy*d      uy*uz*d - ux*s;
        uz*ux*d - uy*s    uz*uy*d + ux*s   c + uz*uz*d
    ]
end

"""
    old_propagate_attitude(w_rad_s, w_dir, tvec) -> Vector{SMatrix}

Copia del algoritmo del antiguo `Propagator.propagate_attitude` (rotación
incremental desde la identidad), para comprobar que el modo `spin` da lo mismo.

# Argumentos
- `w_rad_s::Real`: módulo de ω (rad/s).
- `w_dir::SVector{3,Float64}`: eje de giro (unitario).
- `tvec::AbstractVector`: instantes (s), equiespaciados.

# Devuelve
- `Vector{SMatrix{3,3,Float64,9}}` con DCM cuerpo -> ECI.

# Excepciones
- Ninguna.
"""
function old_propagate_attitude(w_rad_s::Real, w_dir::SVector{3,Float64}, tvec::AbstractVector)
    D = Vector{SMatrix{3,3,Float64,9}}(undef, length(tvec))
    D[1] = SMatrix{3,3,Float64,9}(I)
    ΔD = axis_angle_rotation(w_dir, w_rad_s * (tvec[2] - tvec[1]))
    for k in 2:length(tvec)
        D[k] = ΔD * D[k-1]
    end
    return D
end

# Instantes y trayectoria comunes (0.2 h cada 2 s, como los JSON de ejemplo)
const T    = collect(0.0:2.0:720.0)
const TRAJ = circular_traj(7000.0, 98.0, T)
const SUN  = [normalize(SVector(0.3, -sind(98.0), cosd(98.0))) for _ in T]   # Sol fijo, fuera del plano orbital

# Error máximo entre dos listas de matrices
maxerr(A1, A2) = maximum(maximum(abs.(a - b)) for (a, b) in zip(A1, A2))

@testset "Attitude" begin

    @testset "Convenciones del T6" begin
        # Cuaternión de 90° alrededor de z = matriz elemental C_z(90°) (diap. 19 y 36-37)
        @test A.quat_to_dcm([cosd(45), 0, 0, sind(45)]) ≈ A.elementary_dcm(3, π/2) atol = 1e-12

        # 3-2-1: C_x(φ) C_y(θ) C_z(ψ), matriz explícita de la diap. 26
        ψ, θ, φ = deg2rad.((30.0, 10.0, 5.0))
        C321 = [cos(θ)*cos(ψ)                              cos(θ)*sin(ψ)                              -sin(θ);
                sin(φ)*sin(θ)*cos(ψ)-cos(φ)*sin(ψ)          sin(φ)*sin(θ)*sin(ψ)+cos(φ)*cos(ψ)          sin(φ)*cos(θ);
                cos(φ)*sin(θ)*cos(ψ)+sin(φ)*sin(ψ)          cos(φ)*sin(θ)*sin(ψ)-sin(φ)*cos(ψ)          cos(φ)*cos(θ)]
        @test A.euler_to_dcm([ψ, θ, φ], "321") ≈ C321 atol = 1e-12

        # 3-1-3: C_z(μ) C_x(γ) C_z(ψ), matriz explícita de la diap. 24
        ψ, γ, μr = deg2rad.((40.0, 20.0, 70.0))
        c, s = cos, sin
        C313 = [c(μr)*c(ψ)-s(μr)*s(ψ)*c(γ)    c(μr)*s(ψ)+s(μr)*c(ψ)*c(γ)    s(μr)*s(γ);
                -s(μr)*c(ψ)-c(μr)*s(ψ)*c(γ)   -s(μr)*s(ψ)+c(μr)*c(ψ)*c(γ)   c(μr)*s(γ);
                s(ψ)*s(γ)                     -c(ψ)*s(γ)                    c(γ)]
        @test A.euler_to_dcm([ψ, γ, μr], "313") ≈ C313 atol = 1e-12

        # q y -q son la misma actitud; cuaternión escrito con pocos decimales se acepta
        q = [0.5, 0.5, 0.5, 0.5]
        @test A.quat_to_dcm(q) ≈ A.quat_to_dcm(-q)
        @test A.quat_to_dcm([0.7071, 0, 0, 0.7071]) ≈ A.elementary_dcm(3, π/2) atol = 1e-4
        @test_throws ErrorException A.quat_to_dcm([1, 0, 0, 0.5])         # no unitario
        @test_throws ErrorException A.euler_to_dcm([0, 0, 0], "331")      # dos ejes seguidos iguales
    end

    @testset "spin" begin
        # Compatibilidad: mismo resultado que el antiguo propagate_attitude
        spec = parse_attitude(Dict("type" => "spin", "w_deg_s" => 30, "w_dir" => [1, 0, 0]), "deb")
        @test maxerr(compute_attitude(spec, TRAJ, SUN).DCM_b2i, old_propagate_attitude(deg2rad(30), SVector(1.0, 0, 0), T)) < 1e-10

        # Sin giro: identidad constante (como los JSON antiguos con w_deg_s = 0)
        spec0 = parse_attitude(Dict("type" => "spin", "w_deg_s" => 0), "deb")
        @test all(D -> D ≈ I, compute_attitude(spec0, TRAJ, SUN).DCM_b2i)

        # Ecuación cinemática dC_ib/dt = C_ib [ω_b×] (diferencias centradas con paso h pequeño)
        spec2 = parse_attitude(Dict("type" => "spin", "w_deg_s" => 3, "w_dir" => [1, 2, 2],
                                    "euler_deg" => [10, 20, 30], "sequence" => "321"), "deb")
        h = 1e-3
        Dh = compute_attitude(spec2, circular_traj(7000.0, 98.0, [0.0, 200.0 - h, 200.0, 200.0 + h]), SUN).DCM_b2i
        ωb = deg2rad(3) * SVector(1.0, 2.0, 2.0) / 3
        @test (Dh[4] - Dh[2]) / (2h) ≈ Dh[3] * A.skew(ωb) rtol = 1e-6
        D = compute_attitude(spec2, TRAJ, SUN).DCM_b2i

        # Actitud inicial = transpuesta de la C_bi dada con Euler
        @test D[1] ≈ A.euler_to_dcm(deg2rad.([10, 20, 30]), "321")'

        @test_throws ErrorException parse_attitude(Dict("type" => "spin", "w_deg_s" => 2), "deb")   # falta w_dir
    end

    @testset "random" begin
        b = Dict{String,Any}("type" => "random", "seed" => 7)
        D1 = compute_attitude(parse_attitude(b, "deb"), TRAJ, SUN).DCM_b2i
        D2 = compute_attitude(parse_attitude(b, "deb"), TRAJ, SUN).DCM_b2i
        D3 = compute_attitude(parse_attitude(merge(b, Dict("seed" => 8)), "deb"), TRAJ, SUN).DCM_b2i
        @test maxerr(D1, D2) == 0                       # misma semilla -> mismo resultado
        @test maxerr(D1, D3) > 1e-3                     # otra semilla -> otra actitud

        # |ω| dentro del rango del parámetro maestro, para muchas semillas
        wmin, wmax = deg2rad.(RANDOM_W_RANGE_DEG_S)
        @test all(wmin <= A.sample_random_state(parse_attitude(merge(b, Dict("seed" => s)), "deb")).w_rad_s <= wmax for s in 1:500)

        # El rango ya no se da en el JSON
        @test_throws ErrorException parse_attitude(merge(b, Dict("w_range_deg_s" => [0, 5])), "deb")

        # Orientación uniforme: la media de las DCM sobre todas las rotaciones es la matriz nula
        rng = MersenneTwister(1)
        M = sum(A.quat_to_dcm(A.random_quaternion(rng)) for _ in 1:20000) / 20000
        @test maximum(abs.(M)) < 0.03
    end

    @testset "Propagación del cuaternión" begin
        # Ω² = -|ω|² I₄ (propiedad que da la solución cerrada de la exponencial)
        ω = SVector(0.3, -0.2, 0.5)
        @test A.omega_matrix(ω)^2 ≈ -dot(ω, ω) * I

        # DCM -> cuaternión -> DCM, incluidos giros de ≈ 180° (q_η ≈ 0, donde falla la fórmula directa)
        rng = MersenneTwister(2)
        qs = [A.random_quaternion(rng) for _ in 1:1000]
        @test all(isapprox(A.quat_to_dcm(A.dcm_to_quat(A.quat_to_dcm(q))), A.quat_to_dcm(q); atol = 1e-12) for q in qs)
        q180 = SVector(0.0, 0.0, 1.0, 0.0)                              # 180° alrededor de y
        @test A.dcm_to_quat(A.quat_to_dcm(q180)) ≈ q180

        # Paso a paso = solución analítica de salto único C_ib(t) = C_ib(0)·R(ê, ωt) (ω constante)
        q0, ê, w = A.random_quaternion(rng), A.random_unit_vector(rng), deg2rad(4.0)
        tt = collect(0.0:2.0:86400.0)                                    # un día cada 2 s (43 201 pasos)
        Dstep = A.spin_history(q0, ê, w, tt)
        C0 = A.quat_to_dcm(q0)'
        @test maximum(maximum(abs.(Dstep[k] - C0 * axis_angle_rotation(ê, w * tt[k]))) for k in eachindex(tt)) < 1e-9

        # La norma se mantiene en 1 tras muchos pasos (renormalización)
        q = q0
        for _ in 1:100000
            q = A.quat_step(q, w * ê, 2.0)
        end
        @test abs(norm(q) - 1) < 1e-14
    end

    @testset "fixed" begin
        # ECI: constante
        s = parse_attitude(Dict("type" => "fixed", "frame" => "ECI", "q" => [0.5, 0.5, 0.5, 0.5]), "deb")
        D = compute_attitude(s, TRAJ, SUN).DCM_b2i
        @test all(Dk -> Dk == D[1], D)
        @test D[1] ≈ A.quat_to_dcm([0.5, 0.5, 0.5, 0.5])'

        # LVLH con orientación identidad: ejes cuerpo = (radial, traza, normal)
        s = parse_attitude(Dict("type" => "fixed", "frame" => "LVLH", "q" => [1, 0, 0, 0]), "deb")
        D = compute_attitude(s, TRAJ, SUN).DCM_b2i
        r, v = TRAJ.r[50], TRAJ.v[50]
        @test D[50] * SVector(1, 0, 0) ≈ r / norm(r)
        @test D[50] * SVector(0, 0, 1) ≈ normalize(cross(r, v))
        @test D[50] * SVector(0, 1, 0) ≈ v / norm(v) atol = 1e-12     # órbita circular: traza = velocidad

        @test_throws ErrorException parse_attitude(Dict("type" => "fixed", "frame" => "LVLH"), "deb")          # sin orientación
        @test_throws ErrorException parse_attitude(Dict("type" => "fixed", "frame" => "RTN", "q" => [1, 0, 0, 0]), "deb")
    end

    @testset "pointing" begin
        b = Dict("type" => "pointing",
                 "primary"   => Dict("body_axis" => [0, 0, 1], "direction" => "nadir"),
                 "secondary" => Dict("body_axis" => [1, 0, 0], "direction" => "sun"))
        D = compute_attitude(parse_attitude(b, "sat"), TRAJ, SUN).DCM_b2i
        for k in (1, 120, 361)
            r = TRAJ.r[k]
            @test D[k] * SVector(0, 0, 1) ≈ -r / norm(r)                       # primario exacto
            x_eci = D[k] * SVector(1, 0, 0)
            @test abs(dot(x_eci, -r / norm(r))) < 1e-12                        # x ⟂ z
            @test dot(x_eci, SUN[k]) >= 0                                      # x hacia el lado del Sol
            @test det(D[k]) ≈ 1 && D[k]' * D[k] ≈ I                           # rotación propia
        end

        # Direcciones siempre paralelas o ejes del cuerpo paralelos -> error al leer
        bad = Dict("type" => "pointing",
                   "primary"   => Dict("body_axis" => [0, 0, 1], "direction" => "nadir"),
                   "secondary" => Dict("body_axis" => [1, 0, 0], "direction" => "zenith"))
        @test_throws ErrorException parse_attitude(bad, "sat")
        bad["secondary"] = Dict("body_axis" => [0, 0, -1], "direction" => "sun")
        @test_throws ErrorException parse_attitude(bad, "sat")
    end

    @testset "camera_mount" begin
        # Montaje por defecto: ejes cámara = ejes cuerpo
        @test A.parse_camera_mount(A.DEFAULT_CAMERA_MOUNT, "cam1") ≈ I

        # Cámara mirando a +x cuerpo con +z cuerpo arriba en la imagen
        C_cb = A.parse_camera_mount(Dict("boresight" => [1, 0, 0], "up" => [0, 0, 1]), "cam1")
        @test C_cb * SVector(1, 0, 0) ≈ SVector(0, 0, 1)                  # boresight -> z_c
        @test C_cb * SVector(0, 0, 1) ≈ SVector(0, -1, 0)                 # arriba -> -y_c
        @test det(C_cb) ≈ 1 && C_cb' * C_cb ≈ I                           # rotación propia (T6, diap. 17)

        # "up" no perpendicular: se proyecta sobre el plano del sensor
        C2 = A.parse_camera_mount(Dict("boresight" => [0, 0, 1], "up" => [0, -1, 1]), "cam1")
        @test C2 ≈ I

        # Sin bloque: montaje por defecto (con aviso)
        @test A.camera_mount_block(Dict("SMA" => 7000), "cam1") == A.DEFAULT_CAMERA_MOUNT

        @test_throws ErrorException A.parse_camera_mount(Dict("boresight" => [0, 0, 1], "up" => [0, 0, -2]), "cam1")  # paralelos
        @test_throws ErrorException A.parse_camera_mount(Dict("boresight" => [0, 0, 1], "upp" => [0, 1, 0]), "cam1")  # errata
        @test_throws ErrorException A.parse_camera_mount(Dict("boresight" => [0, 0, 1]), "cam1")                     # falta up
    end

    @testset "tracking" begin
        # Objetivo en una órbita circular algo más alta y desfasada: la línea de visión cambia con el tiempo
        obj = circular_traj(7100.0, 97.0, T .+ 30.0)
        trajs = Dict("deb1" => obj)
        b = Dict("type" => "tracking", "objective" => "deb1",
                 "roll" => Dict("body_axis" => [0, -1, 0], "direction" => [0, 0, 1]))
        spec = parse_attitude(b, "cam1"; role = "Observer")

        for mount in (Dict("boresight" => [0, 0, 1], "up" => [0, -1, 0]),   # cámara por defecto
                      Dict("boresight" => [1, 0, 0], "up" => [0, 0, 1]))    # cámara mirando a +x cuerpo
            C_cb = A.parse_camera_mount(mount, "cam1")
            bound = bind_attitude(spec, (camera_mount = C_cb, objective_names = keys(trajs), trajectories = trajs))
            att = compute_attitude(bound, TRAJ, SUN)
            C_ci = camera_dcm_history(C_cb, att)
            b_body = SVector{3,Float64}(mount["boresight"])
            for k in (1, 180, 361)
                û = normalize(obj.r[k] - TRAJ.r[k])
                @test att.DCM_b2i[k] * b_body ≈ û                              # el boresight apunta al RSO
                @test C_ci[k][3, :] ≈ û                                        # fila 3 de C_ci = boresight en ECI
                @test C_ci[k] * û ≈ SVector(0, 0, 1)                           # el RSO queda en el centro de la imagen
                @test det(C_ci[k]) ≈ 1 && C_ci[k]' * C_ci[k] ≈ I
                # El eje de "roll" queda en el plano formado por la línea de visión y la dirección pedida
                y_eci = att.DCM_b2i[k] * SVector(0, -1, 0)
                @test abs(dot(y_eci, normalize(cross(û, SVector(0, 0, 1))))) < 1e-10
                @test dot(y_eci, SVector(0, 0, 1)) >= 0
            end
        end

        # Errores: sin cámara (objetivo), objetivo desconocido, roll paralelo al boresight, sin bind
        C_cb = A.parse_camera_mount(A.DEFAULT_CAMERA_MOUNT, "cam1")
        @test_throws ErrorException bind_attitude(spec, (camera_mount = nothing, objective_names = nothing, trajectories = nothing))
        @test_throws ErrorException bind_attitude(spec, (camera_mount = C_cb, objective_names = ["deb2"], trajectories = nothing))
        bpar = merge(b, Dict("roll" => Dict("body_axis" => [0, 0, 1], "direction" => [0, 0, 1])))
        @test_throws ErrorException bind_attitude(parse_attitude(bpar, "cam1"; role = "Observer"),
                                                  (camera_mount = C_cb, objective_names = nothing, trajectories = nothing))
        @test_throws ErrorException compute_attitude(spec, TRAJ, SUN)
        @test_throws ErrorException parse_attitude(merge(b, Dict("target" => "deb1")), "cam1"; role = "Observer")   # clave no admitida
    end

    @testset "Orientación de la imagen" begin
        # Cámara por defecto apuntando a ECI +y con el norte celeste (ECI +z) arriba en la imagen
        b = Dict("type" => "pointing",
                 "primary"   => Dict("body_axis" => [0, 0, 1],  "direction" => [0, 1, 0]),
                 "secondary" => Dict("body_axis" => [0, -1, 0], "direction" => [0, 0, 1]))
        att = compute_attitude(parse_attitude(b, "cam1"; role = "Observer"), TRAJ, SUN)
        C_ci = camera_dcm_history(A.parse_camera_mount(A.DEFAULT_CAMERA_MOUNT, "cam1"), att)[1]
        bs = SVector(0.0, 1.0, 0.0)
        east = SVector(-1.0, 0.0, 0.0)          # este en RA = 90°: dirección de RA creciente, ẑ × b̂
        north = SVector(0.0, 0.0, 1.0)

        # Cámara real (dextrógira): con el norte arriba, el este sale a la IZQUIERDA (x_c < 0)
        @test (C_ci * normalize(bs + 0.01 * east))[1] < 0
        @test (C_ci * normalize(bs + 0.01 * north))[2] < 0           # norte hacia arriba (y_c < 0)
        @test det(C_ci) ≈ 1                                         # no es un espejo
    end

    @testset "Lectura del JSON del simulador" begin
        # Los mensajes de error indican si es un objetivo o un observador
        err = try
            parse_member_attitude(Dict("SMA" => 7000), "cam1"; role = "Observer"); ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("Observer 'cam1'", err)

        @test_throws ErrorException parse_member_attitude(Dict("SMA" => 7000), "deb")                     # sin actitud
        @test_throws ErrorException parse_attitude(Dict("type" => "tumbling"), "deb")                    # modo desconocido
        @test_throws ErrorException parse_attitude(Dict("type" => "spin", "w_deg_s" => 0, "w_deg" => 1), "deb")   # errata
    end
end
