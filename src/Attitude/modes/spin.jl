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
# q̇ = ½ Ω(ω_b) q, paso a paso: q_{k+1} = exp(½ Ω Δt_k) q_k, renormalizando en cada
# paso (quat_step). Con ω_b = ω ê constante esta actualización es exacta: equivale a
# C_ib(t) = C_ib(0) · R(ê, ω t). La DCM cuerpo -> ECI de cada instante es C(q)ᵀ.
# No se modela la inercia: el giro solo es dinámicamente consistente si ê es un eje
# principal de inercia (hipótesis aceptada con el tutor: los objetos se ven como puntos).
#
# Funciones: parse_spin, spin_history, compute_attitude(::SpinSpec, ...)
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
- `name::AbstractString`: nombre del objetivo, para los mensajes de error.

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
    who = "Objective '$name' (attitude 'spin')"
    check_keys(b, SPIN_KEYS, who)

    w = getnum(b, "w_deg_s", who)

    # Sin giro el eje no importa; con giro es obligatorio
    axis = if haskey(b, "w_dir")
        getunit(b, "w_dir", who)
    else
        w == 0 || error("$who: 'w_dir' is required when 'w_deg_s' ≠ 0.")
        Vec3(1.0, 0.0, 0.0)
    end

    # Actitud inicial opcional (C_bi) como cuaternión; por defecto, ejes cuerpo = ECI
    C_bi0 = parse_orientation(b, who; required = false)
    q0 = C_bi0 === nothing ? SVector(1.0, 0.0, 0.0, 0.0) : dcm_to_quat(C_bi0)

    return SpinSpec(q0, axis, deg2rad(w))
end

"""
    spin_history(q0, axis_b, w_rad_s, t) -> Vector{SMatrix{3,3,Float64,9}}

Historia de la DCM cuerpo -> ECI de un cuerpo que gira con ω constante en ejes
cuerpo, propagando el cuaternión paso a paso con `quat_step`:
`q_k = exp(½ Ω(ω) (t_k - t_{k-1})) · q_{k-1}` y `DCM_b2i(t_k) = C(q_k)ᵀ`.

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
