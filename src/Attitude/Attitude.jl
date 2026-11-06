module Attitude

# =============================================================================
# Módulo Attitude
#
# Define la actitud de los objetivos (RSOs) y de los satélites observadores a lo
# largo de la simulación. Cada uno lleva en su JSON un bloque "attitude" con un
# modo ("type"); el módulo lo lee, lo valida y calcula la historia de la DCM
# cuerpo -> ECI, que usan la fotometría (objetivos) y la cámara (observadores).
#
#   "attitude": {"type": "spin",     "w_deg_s": 2.0, "w_dir": [0, 0, 1]}
#   "attitude": {"type": "random",   "seed": 42}
#   "attitude": {"type": "fixed",    "frame": "LVLH", "euler_deg": [30, 0, 0], "sequence": "321"}
#   "attitude": {"type": "pointing", "primary":   {"body_axis": [0, 0, 1], "direction": "nadir"},
#                                    "secondary": {"body_axis": [1, 0, 0], "direction": "velocity"}}
#   "attitude": {"type": "tracking", "objective": "deb1",                       (solo observadores)
#                                    "roll": {"body_axis": [0, -1, 0], "direction": [0, 0, 1]}}
#
# Los observadores llevan además un bloque "camera_mount" (Mount.jl) con el montaje de la
# cámara en el satélite (C_cb). La orientación de la cámara es C_ci = C_cb · C_bi.
#
# Scripts que contiene:
#   Attitude.jl     : este archivo; tipo base, registro de modos, lectura del bloque y autocarga
#   Quaternions.jl  : parametrizaciones de actitud (cuaternión, Euler, DCM elementales, sorteo
#                     uniforme) y propagación del cuaternión paso a paso
#   Frames.jl       : sistemas y direcciones de referencia (LVLH, nadir, Sol...) y método TRIAD
#   Mount.jl        : montaje de la cámara en el satélite (bloque "camera_mount", C_cb) y C_ci
#   modes/          : un archivo por modo: fixed.jl, pointing.jl, random.jl, spin.jl, tracking.jl
#
# Convenciones (las del tema T6, Attitude Kinematics):
#   - C_bi: DCM que pasa componentes de ECI a cuerpo (r_b = C_bi r_i). C_ib = C_biᵀ.
#   - Cuaternión escalar primero, q = [q_η, q_1, q_2, q_3], unitario; q y -q son la misma actitud.
#   - C(q) = (q_η² - q_vᵀq_v) I + 2 q_v q_vᵀ - 2 q_η [q_v×]  (= C_bi).
#   - Euler: secuencia "abc" con ángulos [θ1, θ2, θ3] en orden de aplicación:
#     C = C_c(θ3) · C_b(θ2) · C_a(θ1), con C_x, C_y, C_z las matrices elementales del T6.
#   - Lo que se entrega al simulador es AttitudeData.DCM_b2i = C_ib (cuerpo -> ECI), porque la
#     fotometría hace n_eci = DCM_b2i · n_body.
#   - Todos los sistemas de referencia son triedros dextrógiros (T6, diap. 4) y det C = +1 (diap. 17).
#
# Contrato de un modo (para añadir uno nuevo basta con crear un archivo en modes/):
#   1. struct NuevoSpec <: AttitudeSpec ... end                    parámetros ya validados
#   2. parse_nuevo(b::Dict{String,Any}, who) -> NuevoSpec           lee y valida el bloque
#      (who = "Objective 'deb'" u "Observer 'cam1'", para los mensajes de error)
#   3. compute_attitude(spec::NuevoSpec, traj, sun_dir) -> AttitudeData
#   4. register_attitude!("nuevo", parse_nuevo)
#   5. (opcional) bind_attitude(spec::NuevoSpec, ctx) -> NuevoSpec, si el modo necesita
#      datos de fuera de su trayectoria (montaje de la cámara, trayectorias de otros objetos)
#   No hay que tocar nada más: este módulo incluye todos los .jl de modes/.
#
# Dependencias: Types (AttitudeData, TrajectoryData), Parameters (RANDOM_W_RANGE_DEG_S)
# y JSONUtils (lectores de bloques), que deben estar cargados antes.
# =============================================================================

export AttitudeSpec, register_attitude!, parse_attitude, parse_member_attitude
export compute_attitude, attitude_block, bind_attitude
export camera_mount_block, parse_camera_mount, camera_dcm_history

using ..Types: AttitudeData, TrajectoryData
using ..Parameters: RANDOM_W_RANGE_DEG_S         # rango de |ω| del modo "random" (parámetro maestro)
using ..JSONUtils     # lectores de bloques: todict, check_keys, getnum, getvec, getunit, get_or_draw_seed

using LinearAlgebra
using Random
using StaticArrays

# Alias de los tipos que más se repiten
const Vec3 = SVector{3,Float64}
const Mat3 = SMatrix{3,3,Float64,9}

# Tipo base de todos los modos; cada modo define su subtipo con sus parámetros
abstract type AttitudeSpec end

# -----------------------------------------------------------------------------
# Registro de modos
# -----------------------------------------------------------------------------

# Nombre del modo ("spin", "random"...) -> función que lee su bloque
const ATTITUDES = Dict{String,Function}()

# Registra la función que lee el bloque de un modo; devuelve nothing
register_attitude!(name::AbstractString, f::Function) = (ATTITUDES[String(name)] = f; nothing)

# -----------------------------------------------------------------------------
# Lectura del bloque "attitude"
# -----------------------------------------------------------------------------

"""
    parse_attitude(block, name; role="Objective") -> AttitudeSpec

Lee y valida un bloque `attitude` con el modo indicado en su campo `type`.

# Argumentos
- `block::AbstractDict`: bloque `attitude` (de JSON3 o de un `Dict`).
- `name::AbstractString`: nombre del objetivo u observador, para los mensajes de error.
- `role::AbstractString`: `"Objective"` u `"Observer"`; solo cambia el texto de los mensajes.

# Devuelve
- Un subtipo de `AttitudeSpec` (`SpinSpec`, `RandomSpec`, `FixedSpec`, `PointingSpec`, `TrackingSpec`...).

# Excepciones
- `ErrorException`: si el bloque no es un diccionario, no tiene `type` o el modo no existe.
- `ErrorException`: cualquier error de validación del modo elegido.

# Ejemplo
```julia
parse_attitude(Dict("type" => "spin", "w_deg_s" => 2, "w_dir" => [0, 0, 1]), "deb")
parse_attitude(Dict("type" => "fixed", "frame" => "ECI", "q" => [1, 0, 0, 0]), "cam1"; role = "Observer")
```
"""
function parse_attitude(block, name::AbstractString; role::AbstractString = "Objective")
    who = "$role '$name'"                           # p. ej. "Observer 'cam1'"
    b = todict(block)
    b isa Dict{String,Any} || error("$who: 'attitude' must be a block {...} (got $(repr(block))).")
    haskey(b, "type") || error("$who: 'attitude' block without 'type'.")

    # El modo se busca en el registro; si no existe se listan los disponibles
    mode = String(b["type"])
    haskey(ATTITUDES, mode) || error(
        "$who: unknown attitude type '$mode'. Available: " * join(sort!(collect(keys(ATTITUDES))), ", "))

    return ATTITUDES[mode](b, who)
end

"""
    attitude_block(entry, name; role="Objective") -> Dict{String,Any}

Devuelve el bloque `attitude` de un objetivo u observador (del JSON del simulador
o de un grupo de la spec) tras comprobar que existe y es un bloque `{...}`.

Por qué existe: el simulador y el generador aplican exactamente las mismas reglas
y los mismos mensajes de error.

# Argumentos
- `entry::AbstractDict`: entrada del miembro (o grupo de la spec).
- `name::AbstractString`: nombre del miembro, para los mensajes de error.
- `role::AbstractString`: `"Objective"` u `"Observer"`; solo cambia el texto de los mensajes.

# Devuelve
- `Dict{String,Any}`: copia del bloque `attitude` (sin validar sus parámetros).

# Excepciones
- `ErrorException`: si falta el bloque `attitude` o si no es un bloque `{...}`.

# Ejemplo
```julia
attitude_block(Dict("attitude" => Dict("type" => "spin", "w_deg_s" => 0)), "deb")
# Dict("type" => "spin", "w_deg_s" => 0)
```
"""
function attitude_block(entry, name::AbstractString; role::AbstractString = "Objective")
    e = todict(entry)
    haskey(e, "attitude") || error(
        "$role '$name': missing 'attitude' block (e.g. \"attitude\": {\"type\": \"spin\", \"w_deg_s\": 0}).")
    e["attitude"] isa Dict{String,Any} || error("$role '$name': 'attitude' must be a block {...}.")
    return e["attitude"]
end

"""
    parse_member_attitude(entry, name; role="Objective") -> AttitudeSpec

Lee la actitud de un objetivo u observador del JSON del simulador: comprueba el
bloque con `attitude_block` y lo valida con `parse_attitude`.

# Argumentos
- `entry::AbstractDict`: entrada del miembro en la sección `objectives` u `observers`.
- `name::AbstractString`: nombre del miembro, para los mensajes de error.
- `role::AbstractString`: `"Objective"` u `"Observer"`.

# Devuelve
- Un subtipo de `AttitudeSpec`.

# Excepciones
- Las de `attitude_block` y `parse_attitude`.

# Ejemplo
```julia
entry = Dict("SMA" => 7000, "attitude" => Dict("type" => "spin", "w_deg_s" => 0))
parse_member_attitude(entry, "deb")       # SpinSpec sin giro
```
"""
function parse_member_attitude(entry, name::AbstractString; role::AbstractString = "Objective")
    return parse_attitude(attitude_block(entry, name; role = role), name; role = role)
end

# -----------------------------------------------------------------------------
# Cálculo de la actitud
# -----------------------------------------------------------------------------

"""
    bind_attitude(spec, ctx) -> AttitudeSpec

Completa un modo con datos que no están en su propia trayectoria: el montaje de
la cámara o las trayectorias de otros objetos. Se llama antes de
`compute_attitude`. Este método de reserva no hace nada (devuelve `spec`); los
modos que lo necesitan (p. ej. `tracking`) definen el suyo.

Por qué existe: así `compute_attitude(spec, traj, sun_dir)` mantiene la misma
firma en todos los modos y no hay que tocar los modos existentes.

# Argumentos
- `spec::AttitudeSpec`: actitud leída con `parse_attitude`.
- `ctx`: `NamedTuple` con
  - `camera_mount`: C_cb del observador (`SMatrix{3,3}`), o `nothing` si es un objetivo.
  - `objective_names`: nombres de los objetivos del escenario (para validar), o `nothing`.
  - `trajectories`: `Dict` nombre del objetivo -> `TrajectoryData`, o `nothing` (generador).

# Devuelve
- El mismo `spec` o una copia completada.

# Excepciones
- Las del método de cada modo (p. ej. objetivo inexistente en `tracking`).

# Ejemplo
```julia
spec = bind_attitude(spec, (camera_mount = C_cb, objective_names = keys(trajs), trajectories = trajs))
```
"""
function bind_attitude(spec::AttitudeSpec, ctx)
    return spec
end

"""
    compute_attitude(spec, traj, sun_dir) -> AttitudeData

Calcula la historia de la DCM cuerpo -> ECI de un objetivo u observador. Cada
modo añade su propio método; este es el de reserva, que solo da error.

# Argumentos
- `spec::AttitudeSpec`: actitud leída con `parse_attitude` (y completada con `bind_attitude`).
- `traj::TrajectoryData`: trayectoria del satélite (r, v en ECI, km y km/s; t en s).
- `sun_dir::AbstractVector`: dirección unitaria del Sol en ECI en cada instante de `traj.t`.

# Devuelve
- `AttitudeData` con un `DCM_b2i` por instante de `traj.t`.

# Excepciones
- `ErrorException`: si el modo no tiene método `compute_attitude`.
"""
function compute_attitude(spec::AttitudeSpec, traj::TrajectoryData, sun_dir::AbstractVector)
    error("compute_attitude is not implemented for $(typeof(spec)).")
end

# -----------------------------------------------------------------------------
# Utilidades + autocarga de modos
# -----------------------------------------------------------------------------

include(joinpath(@__DIR__, "Quaternions.jl"))
include(joinpath(@__DIR__, "Frames.jl"))
include(joinpath(@__DIR__, "Mount.jl"))          # montaje de la cámara (C_cb) y C_ci

# Se incluyen todos los .jl de modes/, en orden alfabético
let dir = joinpath(@__DIR__, "modes")
    # Si el módulo se precompila, esto permite detectar archivos nuevos; si falla, no importa
    try
        include_dependency(dir)
    catch
    end
    for f in sort!(filter(endswith(".jl"), readdir(dir)))
        include(joinpath(dir, f))
    end
end

end # module Attitude
