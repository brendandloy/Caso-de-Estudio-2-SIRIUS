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
# Modelo: con ω_b = ω ê constante, la ecuación cinemática dC_ib/dt = C_ib [ω_b×]
# tiene solución exacta C_ib(t) = C_ib(0) · R(ê, ω t), con R la rotación de
# Rodrigues (axis_angle_rotation). No es dinámicamente consistente salvo que ê sea
# un eje principal de inercia; para eso hará falta un modo con ecuaciones de Euler.
#
# Funciones: parse_spin, spin_history, compute_attitude(::SpinSpec, ...)
# =============================================================================

"""
    SpinSpec(C_ib0, axis_b, w_rad_s)

Parámetros validados del modo `spin`.

# Campos
- `C_ib0::SMatrix{3,3,Float64,9}`: actitud inicial, cuerpo -> ECI.
- `axis_b::SVector{3,Float64}`: eje de giro unitario en ejes cuerpo.
- `w_rad_s::Float64`: módulo de la velocidad angular (rad/s).
"""
struct SpinSpec <: AttitudeSpec
    C_ib0::Mat3
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

    # Actitud inicial opcional (C_bi); por defecto, ejes cuerpo = ECI
    C_bi0 = parse_orientation(b, who; required = false)
    C_ib0 = C_bi0 === nothing ? Mat3(I) : Mat3(C_bi0')

    return SpinSpec(C_ib0, axis, deg2rad(w))
end

"""
    spin_history(C_ib0, axis_b, w_rad_s, t) -> Vector{SMatrix{3,3,Float64,9}}

Historia de la DCM cuerpo -> ECI de un cuerpo que gira con ω constante:
`C_ib(t) = C_ib(0) · R(axis_b, w_rad_s · (t - t₀))`.

Por qué existe: la comparten los modos `spin` y `random` (este último solo
sortea las condiciones iniciales).

# Argumentos
- `C_ib0::SMatrix{3,3,Float64,9}`: actitud en el primer instante.
- `axis_b::SVector{3,Float64}`: eje de giro unitario en ejes cuerpo.
- `w_rad_s::Real`: módulo de ω (rad/s).
- `t::AbstractVector`: instantes (s).

# Devuelve
- `Vector{SMatrix{3,3,Float64,9}}` con una DCM por instante.

# Excepciones
- Ninguna.

# Ejemplo
```julia
# 90°/s alrededor de z: tras 1 s el eje x del cuerpo apunta a +y (ECI)
C = spin_history(SMatrix{3,3}(1.0I), SVector(0.0, 0.0, 1.0), deg2rad(90), [0.0, 1.0])
C[2] * [1, 0, 0]       # ≈ [0, 1, 0]
```
"""
function spin_history(C_ib0::Mat3, axis_b::Vec3, w_rad_s::Real, t::AbstractVector)
    iszero(w_rad_s) && return fill(C_ib0, length(t))      # sin giro: actitud constante

    t0 = Float64(first(t))
    # Fórmula cerrada en cada instante: no acumula error de redondeo y admite pasos no uniformes
    return [C_ib0 * axis_angle_rotation(axis_b, w_rad_s * (Float64(tk) - t0)) for tk in t]
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
    return AttitudeData(spin_history(spec.C_ib0, spec.axis_b, spec.w_rad_s, traj.t))
end

register_attitude!("spin", parse_spin)
