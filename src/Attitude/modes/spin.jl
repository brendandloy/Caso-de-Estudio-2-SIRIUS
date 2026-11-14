# =============================================================================
# spin.jl
#
# Modo "spin": giro con velocidad angular constante alrededor de un eje fijo
# en el cuerpo. Es el comportamiento de las claves antiguas w_deg_s / w_dir.
#
#   "attitude": {"type": "spin", "w_deg_s": 2.0, "w_dir": [0, 0, 1]}
#   "attitude": {"type": "spin", "w_deg_s": 2.0, "w_dir": [0, 0, 1],
#                "euler_deg": [0, 45, 0], "sequence": "321"}          (actitud inicial)
#
# - "w_deg_s": módulo de ω (°/s). Con 0 el objeto no gira y "w_dir" puede omitirse.
# - "w_dir":   eje de giro en ejes cuerpo (se normaliza).
# - "q" o "euler_deg" + "sequence" (opcional): actitud inicial del cuerpo respecto a ECI,
#   en las convenciones del T6. Sin ella, los ejes cuerpo coinciden con ECI en t = 0.
#
# Modelo: se propaga el cuaternión con la ecuación cinemática del T6 (diap. 49),
# q̇ = ½ Ω(ω_b) q, con RK4 entre cada par de instantes de salida, en subpasos de como
# mucho ATTITUDE_MAX_STEP_DEG y renormalizando tras cada uno (quat_step). Con ω_b = ω ê
# constante la solución exacta es C_ib(t) = C_ib(0) · R(ê, ω t); el RK4 la reproduce con
# un error ≈ 3e-10 rad por vuelta. La DCM cuerpo -> ECI de cada instante es C(q)ᵀ.
# No se modela la inercia: el giro solo es dinámicamente consistente si ê es un eje
# principal de inercia (hipótesis aceptada con el tutor: los objetos se ven como puntos).
#
# Funciones: parse_spin, initial_quaternion, spin_history, compute_attitude(::SpinSpec, ...)
# =============================================================================

"""
    SpinSpec(q0, axis_b, w_rad_s)

Parámetros validados del modo `spin`.

# Campos
- `q0::SVector{4,Float64}`: actitud inicial (cuaternión de C_bi, escalar primero).
- `axis_b::SVector{3,Float64}`: eje de giro unitario en ejes cuerpo.
- `w_rad_s::Float64`: módulo de la velocidad angular (rad/s).
"""
struct SpinSpec <: AttitudeSpec
    q0::SVector{4,Float64}
    axis_b::Vec3
    w_rad_s::Float64
end

# Claves admitidas en el bloque "spin"
const SPIN_KEYS = ("type", "w_deg_s", "w_dir", "q", "euler_deg", "sequence")

"""
    parse_spin(b, name) -> SpinSpec

Lee y valida un bloque `attitude` de tipo `spin`.

# Argumentos
- `b::Dict{String,Any}`: bloque `attitude`.
- `name::AbstractString`: miembro que lo usa (p. ej. `"Objective 'deb'"`), para los mensajes de error.

# Devuelve
- `SpinSpec`.

# Excepciones
- `ErrorException`: si hay claves no admitidas, falta `w_deg_s`, falta `w_dir` con
  `w_deg_s ≠ 0`, `w_dir` tiene norma cero o la actitud inicial no es válida.

# Ejemplo
```julia
parse_spin(Dict{String,Any}("type" => "spin", "w_deg_s" => 2, "w_dir" => [0, 0, 1]), "deb")
```
"""
function parse_spin(b::Dict{String,Any}, name::AbstractString)
    who = "$name (attitude 'spin')"
    check_keys(b, SPIN_KEYS; who = who)

    w = getnum(b, "w_deg_s", who)

    # Sin giro el eje no importa; con giro es obligatorio
    axis = if haskey(b, "w_dir")
        getunit(b, "w_dir", who)
    else
        w == 0 || error("$who: 'w_dir' is required when 'w_deg_s' ≠ 0.")
        Vec3(1.0, 0.0, 0.0)
    end

    # Actitud inicial opcional (C_bi); parse_orientation valida el bloque
    C_bi0 = parse_orientation(b, who; required = false)
    q0 = initial_quaternion(b, C_bi0, who)

    return SpinSpec(q0, axis, deg2rad(w))
end

"""
    initial_quaternion(b, C_bi0, who) -> SVector{4,Float64}

Cuaternión inicial del modo `spin`: el `q` del bloque tal cual (normalizado), o el de la
DCM con la fórmula del T6 (`dcm_to_quat`) si la actitud se dio con ángulos de Euler. Sin
actitud inicial, la identidad (ejes cuerpo = ECI).

Por qué existe: el `q` del usuario no se convierte a DCM y de vuelta, así un giro de 180°
escrito como cuaternión no pasa por la fórmula del T6, que divide por q_η ≈ 0.

# Argumentos
- `b::Dict{String,Any}`: bloque `attitude` (ya validado por `parse_orientation`).
- `C_bi0`: DCM C_bi inicial, o `nothing` si no se dio.
- `who::AbstractString`: miembro y modo, para los mensajes de error.

# Devuelve
- `SVector{4,Float64}` unitario (C_bi, escalar primero).

# Excepciones
- `ErrorException`: si los ángulos de Euler dan un giro a menos de ≈ 0.1° de 180°.
"""
function initial_quaternion(b::Dict{String,Any}, C_bi0, who::AbstractString)
    C_bi0 === nothing && return SVector(1.0, 0.0, 0.0, 0.0)
    if haskey(b, "q")
        q = SVector{4,Float64}(getvec(b, "q", 4, who))
        return q / norm(q)
    end
    try
        return dcm_to_quat(C_bi0)
    catch err
        err isa ErrorException ? error("$who: " * err.msg) : rethrow()
    end
end

"""
    spin_history(q0, axis_b, w_rad_s, t) -> Vector{SMatrix{3,3,Float64,9}}

Historia de la DCM cuerpo -> ECI de un cuerpo que gira con ω constante en ejes
cuerpo, propagando el cuaternión entre instantes con `quat_step` (RK4 con subpasos)
y `DCM_b2i(t_k) = C(q_k)ᵀ`.

Por qué existe: la comparten los modos `spin` y `random` (este último solo
sortea las condiciones iniciales).

# Argumentos
- `q0::SVector{4,Float64}`: actitud en el primer instante (cuaternión de C_bi).
- `axis_b::SVector{3,Float64}`: eje de giro unitario en ejes cuerpo.
- `w_rad_s::Real`: módulo de ω (rad/s).
- `t::AbstractVector`: instantes (s); admite pasos no uniformes.

# Devuelve
- `Vector{SMatrix{3,3,Float64,9}}` con una DCM cuerpo -> ECI por instante.

# Excepciones
- Ninguna.

# Ejemplo
```julia
# 90°/s alrededor de z: tras 1 s el eje x del cuerpo apunta a +y (ECI)
C = spin_history(SVector(1.0, 0, 0, 0), SVector(0.0, 0.0, 1.0), deg2rad(90), [0.0, 1.0])
C[2] * [1, 0, 0]       # ≈ [0, 1, 0]
```
"""
function spin_history(q0::SVector{4,Float64}, axis_b::Vec3, w_rad_s::Real, t::AbstractVector)
    ω_b = w_rad_s * axis_b                                 # velocidad angular en ejes cuerpo (rad/s)
    DCM_b2i = Vector{Mat3}(undef, length(t))

    q = q0
    DCM_b2i[1] = Mat3(quat_to_dcm(q)')                     # C_ib = C_biᵀ
    for k in 2:length(t)
        q = quat_step(q, ω_b, Float64(t[k]) - Float64(t[k-1]))   # un paso de la cinemática
        DCM_b2i[k] = Mat3(quat_to_dcm(q)')
    end
    return DCM_b2i
end

"""
    compute_attitude(spec::SpinSpec, traj, sun_dir) -> AttitudeData

Actitud del modo `spin` en los instantes de la trayectoria. No usa ni la
órbita ni el Sol.

# Argumentos
- `spec::SpinSpec`: parámetros del modo.
- `traj::TrajectoryData`: trayectoria del objetivo (solo se usan los tiempos).
- `sun_dir::AbstractVector`: dirección del Sol (no se usa).

# Devuelve
- `AttitudeData`.

# Excepciones
- Ninguna.
"""
function compute_attitude(spec::SpinSpec, traj::TrajectoryData, sun_dir::AbstractVector)
    return AttitudeData(spin_history(spec.q0, spec.axis_b, spec.w_rad_s, traj.t))
end

register_attitude!("spin", parse_spin)
