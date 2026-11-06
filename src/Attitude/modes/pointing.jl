# =============================================================================
# pointing.jl
#
# Modo "pointing": el satélite apunta un eje del cuerpo a una dirección (nadir,
# Sol, velocidad...) y fija el giro alrededor de ese eje con una segunda
# dirección. Típico de satélites operativos.
#
#   "attitude": {"type": "pointing",
#                "primary":   {"body_axis": [0, 0, 1], "direction": "nadir"},
#                "secondary": {"body_axis": [1, 0, 0], "direction": "velocity"}}
#
# - "primary": el eje "body_axis" (ejes cuerpo) apunta exactamente a "direction".
# - "secondary": "body_axis" se acerca todo lo posible a su "direction" (solo fija el
#   giro alrededor del eje primario). Los dos ejes del cuerpo no pueden ser paralelos.
# - "direction": "nadir", "zenith", "velocity", "orbit_normal", "sun" (con "-" delante
#   para invertir) o un vector ECI [x, y, z]. Ver Frames.jl.
#
# Modelo: en cada instante se evalúan las dos direcciones con r, v y el Sol, y la
# actitud se obtiene con el método TRIAD (Frames.jl). Si en algún instante las dos
# direcciones son paralelas (p. ej. Sol alineado con el nadir) la actitud queda
# indefinida y se da error.
#
# Funciones: parse_pointing, compute_attitude(::PointingSpec, ...)
# =============================================================================

"""
    PointingSpec(b1, b2, d1, d2)

Parámetros validados del modo `pointing`.

# Campos
- `b1::SVector{3,Float64}`, `b2::SVector{3,Float64}`: ejes primario y secundario (ejes cuerpo, unitarios).
- `d1`, `d2`: direcciones primaria y secundaria (`NamedTuple` de `parse_direction`).
"""
struct PointingSpec{T1,T2} <: AttitudeSpec
    b1::Vec3
    b2::Vec3
    d1::T1
    d2::T2
end

# Claves admitidas en el bloque "pointing" y en sus sub-bloques
const POINTING_KEYS = ("type", "primary", "secondary")
const POINTING_AXIS_KEYS = ("body_axis", "direction")

"""
    parse_pointing(b, name) -> PointingSpec

Lee y valida un bloque `attitude` de tipo `pointing`.

# Argumentos
- `b::Dict{String,Any}`: bloque `attitude`.
- `name::AbstractString`: miembro que lo usa (p. ej. `"Objective 'deb'"`), para los mensajes de error.

# Devuelve
- `PointingSpec`.

# Excepciones
- `ErrorException`: si hay claves no admitidas, falta `primary`/`secondary` o alguno de
  sus campos, los ejes del cuerpo son paralelos o las dos direcciones lo son siempre
  (p. ej. "nadir" y "zenith").

# Ejemplo
```julia
b = Dict{String,Any}("type" => "pointing",
    "primary"   => Dict{String,Any}("body_axis" => [0, 0, 1], "direction" => "nadir"),
    "secondary" => Dict{String,Any}("body_axis" => [1, 0, 0], "direction" => "velocity"))
parse_pointing(b, "sat")
```
"""
function parse_pointing(b::Dict{String,Any}, name::AbstractString)
    who = "$name (attitude 'pointing')"
    check_keys(b, POINTING_KEYS; who = who)

    # Lee un sub-bloque {"body_axis": [...], "direction": ...}
    function read_axis(key)
        haskey(b, key) || error("$who: missing '$key' block {\"body_axis\": [x, y, z], \"direction\": ...}.")
        sub = b[key]
        sub isa AbstractDict || error("$who: '$key' must be a block {\"body_axis\": [x, y, z], \"direction\": ...}.")
        w = "$who, '$key'"
        check_keys(sub, POINTING_AXIS_KEYS; who = w)
        haskey(sub, "direction") || error("$w: missing 'direction'.")
        return getunit(sub, "body_axis", w), parse_direction(sub["direction"], w)
    end

    b1, d1 = read_axis("primary")
    b2, d2 = read_axis("secondary")

    # Los ejes del cuerpo deben formar un plano
    norm(cross(b1, b2)) > 1e-8 || error("$who: primary and secondary 'body_axis' are parallel.")

    # Direcciones que son paralelas en todo instante (misma dirección con o sin signo)
    always_parallel = d1.kind === d2.kind && d1.kind !== :inertial ||
                      (d1.kind === :nadir && d2.kind === :zenith) || (d1.kind === :zenith && d2.kind === :nadir) ||
                      (d1.kind === :inertial && d2.kind === :inertial && norm(cross(d1.vec, d2.vec)) <= 1e-8)
    always_parallel && error("$who: primary and secondary directions are parallel; choose two different directions.")

    return PointingSpec(b1, b2, d1, d2)
end

"""
    compute_attitude(spec::PointingSpec, traj, sun_dir) -> AttitudeData

Actitud del modo `pointing`: en cada instante evalúa las dos direcciones y
aplica el método TRIAD.

# Argumentos
- `spec::PointingSpec`: parámetros del modo.
- `traj::TrajectoryData`: trayectoria del objetivo (r, v en ECI).
- `sun_dir::AbstractVector`: dirección unitaria del Sol en ECI, una por instante de `traj.t`.

# Devuelve
- `AttitudeData`.

# Excepciones
- `ErrorException`: si `sun_dir` no tiene un valor por instante, o si en algún
  instante las dos direcciones son paralelas (el mensaje indica el tiempo).
"""
function compute_attitude(spec::PointingSpec, traj::TrajectoryData, sun_dir::AbstractVector)
    N = length(traj.t)
    length(sun_dir) == N || error("pointing: sun_dir has $(length(sun_dir)) samples but the trajectory has $N.")

    DCM_b2i = Vector{Mat3}(undef, N)
    for k in 1:N
        r, v, s = traj.r[k], traj.v[k], sun_dir[k]
        r1 = direction_vector(spec.d1, r, v, s)
        r2 = direction_vector(spec.d2, r, v, s)
        # Se añade el instante al error de TRIAD para saber cuándo ocurre
        DCM_b2i[k] = try
            triad(spec.b1, spec.b2, r1, r2)
        catch err
            error("pointing: at t = $(traj.t[k]) s: " * sprint(showerror, err))
        end
    end
    return AttitudeData(DCM_b2i)
end

register_attitude!("pointing", parse_pointing)
