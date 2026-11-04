# =============================================================================
# tracking.jl
#
# Modo "tracking" (solo observadores): el satélite gira para que el boresight de
# su cámara apunte en todo momento a un RSO concreto (un objetivo del escenario).
#
#   "attitude": {"type": "tracking", "objective": "deb1",
#                "roll": {"body_axis": [0, -1, 0], "direction": [0, 0, 1]}}
#
# - "objective": nombre del objetivo tal como aparece en el JSON completo (p. ej. "deb1",
#   "deb_random2"). El eje que apunta al objetivo no se escribe: es el "boresight" del
#   bloque "camera_mount" del observador (Mount.jl).
# - "roll": fija el giro del satélite alrededor del eje que apunta al RSO. Su "body_axis"
#   (ejes cuerpo) se acerca todo lo posible a su "direction" ("nadir", "sun", "orbit_normal",
#   vector ECI... ver Frames.jl). No puede ser paralelo al boresight.
#
# Modelo, en cada instante t_k:
#   û(t_k) = (r_obj - r_obs) / |r_obj - r_obs|      línea de visión observador -> RSO (ECI)
#   C_ib(t_k) = TRIAD(b̂_cam, roll.body_axis; û, roll.direction)
# donde b̂_cam es el boresight en ejes cuerpo (fila 3 de C_cb). Así C_ib·b̂_cam = û exactamente.
# Supuestos: posiciones reales (no las de navegación nav_), sin tiempo de luz ni aberración;
# el satélite apunta al RSO aunque esté eclipsado u oculto por la Tierra (eso lo decide la
# fotometría). Si en algún instante û es paralela a roll.direction, la actitud queda
# indefinida y se da error.
#
# El modo necesita datos de fuera de la trayectoria del observador (montaje de la cámara y
# trayectoria del objetivo): se los da bind_attitude antes de compute_attitude.
#
# Funciones: parse_tracking, bind_attitude(::TrackingSpec, ...), compute_attitude(::TrackingSpec, ...)
# =============================================================================

"""
    TrackingSpec(who, objective, roll_axis, roll_dir, boresight_b, objective_traj)

Parámetros del modo `tracking`. Los dos últimos campos valen `nothing` hasta que
`bind_attitude` los completa.

# Campos
- `who::String`: miembro que lo usa (p. ej. `"Observer 'cam1'"`), para los mensajes de error.
- `objective::String`: nombre del objetivo que se sigue.
- `roll_axis::SVector{3,Float64}`: eje cuerpo que fija el giro (unitario).
- `roll_dir`: dirección a la que se acerca `roll_axis` (`NamedTuple` de `parse_direction`).
- `boresight_b::Union{Nothing,SVector{3,Float64}}`: boresight de la cámara en ejes cuerpo.
- `objective_traj::Union{Nothing,TrajectoryData}`: trayectoria del objetivo.
"""
struct TrackingSpec{T} <: AttitudeSpec
    who::String
    objective::String
    roll_axis::Vec3
    roll_dir::T
    boresight_b::Union{Nothing,Vec3}
    objective_traj::Union{Nothing,TrajectoryData}
end

# Claves admitidas en el bloque "tracking" y en su sub-bloque "roll"
const TRACKING_KEYS = ("type", "objective", "roll")
const TRACKING_ROLL_KEYS = ("body_axis", "direction")

"""
    parse_tracking(b, name) -> TrackingSpec

Lee y valida un bloque `attitude` de tipo `tracking`. La comprobación del montaje
de la cámara y del nombre del objetivo la hace `bind_attitude`.

# Argumentos
- `b::Dict{String,Any}`: bloque `attitude`.
- `name::AbstractString`: miembro que lo usa (p. ej. `"Observer 'cam1'"`), para los mensajes de error.

# Devuelve
- `TrackingSpec` sin completar (sin boresight ni trayectoria del objetivo).

# Excepciones
- `ErrorException`: si hay claves no admitidas, falta `objective` o no es un texto, falta
  `roll` o alguno de sus campos, o la dirección no es válida.

# Ejemplo
```julia
b = Dict{String,Any}("type" => "tracking", "objective" => "deb1",
    "roll" => Dict{String,Any}("body_axis" => [0, -1, 0], "direction" => [0, 0, 1]))
parse_tracking(b, "Observer 'cam1'")
```
"""
function parse_tracking(b::Dict{String,Any}, name::AbstractString)
    who = "$name (attitude 'tracking')"
    check_keys(b, TRACKING_KEYS, who)

    # Objetivo que se sigue: nombre del miembro en la sección "objectives"
    haskey(b, "objective") || error("$who: missing 'objective' (name of the objective to track).")
    obj = b["objective"]
    (obj isa AbstractString && !isempty(obj)) || error("$who: 'objective' must be a name such as \"deb1\" (got $(repr(obj))).")

    # Sub-bloque "roll": eje del cuerpo y dirección a la que se acerca
    haskey(b, "roll") || error("$who: missing 'roll' block {\"body_axis\": [x, y, z], \"direction\": ...}.")
    roll = b["roll"]
    roll isa AbstractDict || error("$who: 'roll' must be a block {\"body_axis\": [x, y, z], \"direction\": ...}.")
    w = "$who, 'roll'"
    check_keys(roll, TRACKING_ROLL_KEYS, w)
    haskey(roll, "direction") || error("$w: missing 'direction'.")

    return TrackingSpec(String(name), String(obj), getunit(roll, "body_axis", w),
                        parse_direction(roll["direction"], w), nothing, nothing)
end

"""
    bind_attitude(spec::TrackingSpec, ctx) -> TrackingSpec

Completa el modo `tracking` con el boresight de la cámara (de `ctx.camera_mount`)
y, si se da, con la trayectoria del objetivo (de `ctx.trajectories`). Comprueba
también que el objetivo existe.

# Argumentos
- `spec::TrackingSpec`: modo leído con `parse_tracking`.
- `ctx`: `NamedTuple` con `camera_mount` (C_cb o `nothing`), `objective_names`
  (colección de nombres o `nothing`) y `trajectories` (`Dict` nombre -> `TrajectoryData` o `nothing`).

# Devuelve
- `TrackingSpec` con `boresight_b` y, si había trayectorias, `objective_traj`.

# Excepciones
- `ErrorException`: si no hay cámara (el modo se usa en un objetivo), si `roll.body_axis`
  es paralelo al boresight o si el objetivo no existe (el mensaje lista los disponibles).

# Ejemplo
```julia
spec = bind_attitude(spec, (camera_mount = C_cb, objective_names = keys(trajs), trajectories = trajs))
```
"""
function bind_attitude(spec::TrackingSpec, ctx)
    who = "$(spec.who) (attitude 'tracking')"

    # Sin cámara no hay boresight que apuntar: el modo es solo para observadores
    ctx.camera_mount === nothing && error("$who: 'tracking' is only available for observers (it points the camera).")
    C_cb = Mat3(ctx.camera_mount)
    b_cam = Vec3(C_cb[3, 1], C_cb[3, 2], C_cb[3, 3])            # boresight en ejes cuerpo (fila 3 de C_cb)
    norm(cross(b_cam, spec.roll_axis)) > 1e-8 ||
        error("$who: 'roll.body_axis' is parallel to the camera boresight; choose another body axis.")

    # El objetivo debe existir en el escenario
    names = ctx.objective_names
    if names !== nothing
        spec.objective in names || error("$who: unknown objective '$(spec.objective)'. Available: " *
                                         join(sort!(String.(collect(names))), ", "))
    end

    trajs = ctx.trajectories
    traj = trajs === nothing ? nothing : trajs[spec.objective]
    return TrackingSpec(spec.who, spec.objective, spec.roll_axis, spec.roll_dir, b_cam, traj)
end

"""
    compute_attitude(spec::TrackingSpec, traj, sun_dir) -> AttitudeData

Actitud del modo `tracking`: en cada instante calcula la línea de visión al
objetivo y aplica el método TRIAD con el boresight de la cámara como eje primario.

# Argumentos
- `spec::TrackingSpec`: modo completado con `bind_attitude`.
- `traj::TrajectoryData`: trayectoria del observador (r, v en ECI).
- `sun_dir::AbstractVector`: dirección unitaria del Sol en ECI, una por instante de `traj.t`.

# Devuelve
- `AttitudeData`.

# Excepciones
- `ErrorException`: si `spec` no se ha completado con `bind_attitude`, si las trayectorias
  o `sun_dir` no tienen el mismo número de instantes, si observador y objetivo coinciden
  o si en algún instante la línea de visión es paralela a `roll.direction` (el mensaje indica el tiempo).
"""
function compute_attitude(spec::TrackingSpec, traj::TrajectoryData, sun_dir::AbstractVector)
    who = "$(spec.who) (attitude 'tracking')"
    (spec.boresight_b === nothing || spec.objective_traj === nothing) &&
        error("$who: call bind_attitude with the camera mount and the trajectories before compute_attitude.")

    N = length(traj.t)
    obj = spec.objective_traj
    length(obj.t) == N || error("$who: objective '$(spec.objective)' has $(length(obj.t)) samples but the observer has $N.")
    length(sun_dir) == N || error("$who: sun_dir has $(length(sun_dir)) samples but the trajectory has $N.")

    DCM_b2i = Vector{Mat3}(undef, N)
    for k in 1:N
        r, v, s = traj.r[k], traj.v[k], sun_dir[k]
        rel = obj.r[k] - r                                     # observador -> objetivo (km, ECI)
        d = norm(rel)
        d > 0 || error("$who: observer and objective '$(spec.objective)' coincide at t = $(traj.t[k]) s.")
        û = rel / d                                            # línea de visión
        r2 = direction_vector(spec.roll_dir, r, v, s)          # dirección de "roll" en ECI
        # Se añade el instante al error de TRIAD para saber cuándo ocurre
        DCM_b2i[k] = try
            triad(spec.boresight_b, spec.roll_axis, û, r2)
        catch err
            error("$who: at t = $(traj.t[k]) s the line of sight is parallel to 'roll.direction'.")
        end
    end
    return AttitudeData(DCM_b2i)
end

register_attitude!("tracking", parse_tracking)
