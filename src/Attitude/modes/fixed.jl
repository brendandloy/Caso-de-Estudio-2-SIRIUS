# =============================================================================
# fixed.jl
#
# Modo "fixed": orientación constante respecto a ECI o respecto a LVLH.
#
#   "attitude": {"type": "fixed", "frame": "ECI",  "q": [1, 0, 0, 0]}
#   "attitude": {"type": "fixed", "frame": "LVLH", "euler_deg": [0, 30, 0], "sequence": "321"}
#
# - "frame": "ECI" (no gira respecto a las estrellas) o "LVLH" (mantiene la
#   orientación respecto a la órbita; ver Frames.jl para los ejes).
# - "q" o "euler_deg" + "sequence": orientación del cuerpo respecto a "frame", en las
#   convenciones del T6 (C_bf, de "frame" a cuerpo). Es obligatoria.
#
# Modelo:
#   ECI:  C_ib = C_bfᵀ                      (constante)
#   LVLH: C_ib(t) = C_io(t) · C_bfᵀ          (C_io: LVLH -> ECI en cada instante)
#
# Funciones: parse_fixed, compute_attitude(::FixedSpec, ...)
# =============================================================================

"""
    FixedSpec(frame, C_fb)

Parámetros validados del modo `fixed`.

# Campos
- `frame::Symbol`: `:ECI` o `:LVLH`.
- `C_fb::SMatrix{3,3,Float64,9}`: DCM cuerpo -> `frame` (= C_bfᵀ).
"""
struct FixedSpec <: AttitudeSpec
    frame::Symbol
    C_fb::Mat3
end

# Claves admitidas en el bloque "fixed"
const FIXED_KEYS = ("type", "frame", "q", "euler_deg", "sequence")

# Sistemas de referencia admitidos en "frame"
const FIXED_FRAMES = ("ECI", "LVLH")

"""
    parse_fixed(b, name) -> FixedSpec

Lee y valida un bloque `attitude` de tipo `fixed`.

# Argumentos
- `b::Dict{String,Any}`: bloque `attitude`.
- `name::AbstractString`: nombre del objetivo, para los mensajes de error.

# Devuelve
- `FixedSpec`.

# Excepciones
- `ErrorException`: si hay claves no admitidas, falta `frame` o no es ECI/LVLH,
  o falta la orientación o no es válida.

# Ejemplo
```julia
parse_fixed(Dict{String,Any}("type" => "fixed", "frame" => "LVLH",
                             "euler_deg" => [0, 30, 0], "sequence" => "321"), "deb")
```
"""
function parse_fixed(b::Dict{String,Any}, name::AbstractString)
    who = "Objective '$name' (attitude 'fixed')"
    check_keys(b, FIXED_KEYS, who)

    haskey(b, "frame") || error("$who: missing 'frame' (\"ECI\" or \"LVLH\").")
    frame = b["frame"]
    (frame isa AbstractString && frame in FIXED_FRAMES) ||
        error("$who: 'frame' must be \"ECI\" or \"LVLH\" (got $(repr(frame))).")

    C_bf = parse_orientation(b, who; required = true)
    return FixedSpec(Symbol(frame), Mat3(C_bf'))
end

"""
    compute_attitude(spec::FixedSpec, traj, sun_dir) -> AttitudeData

Actitud del modo `fixed` en los instantes de la trayectoria.

# Argumentos
- `spec::FixedSpec`: parámetros del modo.
- `traj::TrajectoryData`: trayectoria del objetivo (r y v se usan con LVLH).
- `sun_dir::AbstractVector`: dirección del Sol (no se usa).

# Devuelve
- `AttitudeData`.

# Excepciones
- `ErrorException`: con LVLH, si en algún instante r y v son paralelos.
"""
function compute_attitude(spec::FixedSpec, traj::TrajectoryData, sun_dir::AbstractVector)
    # ECI: la misma DCM en todos los instantes
    spec.frame === :ECI && return AttitudeData(fill(spec.C_fb, length(traj.t)))

    # LVLH: se compone con el giro del sistema orbital en cada instante
    return AttitudeData([lvlh_to_eci(traj.r[k], traj.v[k]) * spec.C_fb for k in eachindex(traj.t)])
end

register_attitude!("fixed", parse_fixed)
