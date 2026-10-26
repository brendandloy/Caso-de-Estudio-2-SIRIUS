module Attitude

# =============================================================================
# Módulo Attitude
#
# Define la actitud de los objetivos (RSOs) a lo largo de la simulación. Cada
# objetivo lleva en su JSON un bloque "attitude" con un modo ("type"); el
# módulo lo lee, lo valida y calcula la historia de la DCM cuerpo -> ECI que
# usa la fotometría.
#
#   "attitude": {"type": "spin",     "w_deg_s": 2.0, "w_dir": [0, 0, 1]}
#   "attitude": {"type": "random",   "w_range_deg_s": [0, 5], "seed": 42}
#   "attitude": {"type": "fixed",    "frame": "LVLH", "euler_deg": [30, 0, 0], "sequence": "321"}
#   "attitude": {"type": "pointing", "primary":   {"body_axis": [0, 0, 1], "target": "nadir"},
#                                    "secondary": {"body_axis": [1, 0, 0], "target": "velocity"}}
#
# Scripts que contiene:
#   Attitude.jl     : este archivo; tipo base, registro de modos, lectura del bloque y autocarga
#   Quaternions.jl  : parametrizaciones de actitud (cuaternión, Euler, DCM elementales, sorteo uniforme)
#   Frames.jl       : sistemas y direcciones de referencia (LVLH, nadir, Sol...) y método TRIAD
#   modes/          : un archivo por modo: fixed.jl, pointing.jl, random.jl, spin.jl
#
# Convenciones (las del tema T6, Attitude Kinematics):
#   - C_bi: DCM que pasa componentes de ECI a cuerpo (r_b = C_bi r_i). C_ib = C_biᵀ.
#   - Cuaternión escalar primero, q = [q_η, q_1, q_2, q_3], unitario; q y -q son la misma actitud.
#   - C(q) = (q_η² - q_vᵀq_v) I + 2 q_v q_vᵀ - 2 q_η [q_v×]  (= C_bi).
#   - Euler: secuencia "abc" con ángulos [θ1, θ2, θ3] en orden de aplicación:
#     C = C_c(θ3) · C_b(θ2) · C_a(θ1), con C_x, C_y, C_z las matrices elementales del T6.
#   - Lo que se entrega al simulador es AttitudeData.DCM_b2i = C_ib (cuerpo -> ECI), porque la
#     fotometría hace n_eci = DCM_b2i · n_body.
#
# Contrato de un modo (para añadir uno nuevo basta con crear un archivo en modes/):
#   1. struct NuevoSpec <: AttitudeSpec ... end                    parámetros ya validados
#   2. parse_nuevo(b::Dict{String,Any}, name) -> NuevoSpec          lee y valida el bloque
#   3. compute_attitude(spec::NuevoSpec, traj, sun_dir) -> AttitudeData
#   4. register_attitude!("nuevo", parse_nuevo)
#   No hay que tocar nada más: este módulo incluye todos los .jl de modes/.
#
# Dependencias: Types (AttitudeData, TrajectoryData) y Transformations
# (axis_angle_rotation), que deben estar cargados antes.
# =============================================================================

export AttitudeSpec, register_attitude!, parse_attitude, parse_member_attitude
export compute_attitude, attitude_block

using ..Types: AttitudeData, TrajectoryData
using ..Transformations: axis_angle_rotation    # rotación de Rodrigues, la misma que usaba propagate_attitude

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
# Helpers de lectura (comunes a todos los modos)
# -----------------------------------------------------------------------------

"""
    todict(x)

Convierte recursivamente un objeto JSON en tipos base de Julia: diccionarios a
`Dict{String,Any}` y vectores a `Vector{Any}`.

Por qué existe: el simulador lee el JSON con JSON3 (claves `Symbol`) y el
generador con `Dict{String,Any}`. Así los modos leen siempre el mismo tipo.

# Argumentos
- `x`: diccionario, vector o escalar.

# Devuelve
- La misma estructura con claves `String`.

# Excepciones
- Ninguna.

# Ejemplo
```julia
todict(JSON3.read("{\\"type\\": \\"spin\\"}"))     # Dict{String,Any}("type" => "spin")
```
"""
function todict(x)
    x isa AbstractDict   && return Dict{String,Any}(String(k) => todict(v) for (k, v) in x)
    x isa AbstractVector && return Any[todict(v) for v in x]
    return x
end

"""
    check_keys(b, allowed, who)

Comprueba que un bloque no tiene claves fuera de las admitidas.

Por qué existe: una errata (p. ej. `w_deg` en vez de `w_deg_s`) se ignoraría en
silencio y la simulación usaría otra actitud sin avisar.

# Argumentos
- `b::AbstractDict`: bloque leído del JSON.
- `allowed`: colección de claves admitidas.
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si aparece alguna clave no admitida.

# Ejemplo
```julia
check_keys(Dict("type" => "spin", "w_deg" => 1), ("type", "w_deg_s"), "deb")
# error: deb: unknown key(s) w_deg. Allowed: type, w_deg_s
```
"""
function check_keys(b::AbstractDict, allowed, who::AbstractString)
    bad = sort!([k for k in keys(b) if !(k in allowed)])
    isempty(bad) || error("$who: unknown key(s) $(join(bad, ", ")). Allowed: $(join(allowed, ", "))")
    return nothing
end

"""
    getnum(b, key, who) -> Float64

Lee un número obligatorio de un bloque.

# Argumentos
- `b::AbstractDict`: bloque leído del JSON.
- `key::AbstractString`: nombre del parámetro.
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.

# Devuelve
- `Float64`: el valor leído.

# Excepciones
- `ErrorException`: si falta, no es un número o no es finito.
"""
function getnum(b::AbstractDict, key::AbstractString, who::AbstractString)
    haskey(b, key) || error("$who: missing '$key'.")
    x = b[key]
    (x isa Real && !(x isa Bool)) || error("$who: '$key' must be a number (got $(repr(x))).")
    isfinite(x) || error("$who: '$key' must be finite.")
    return Float64(x)
end

"""
    getvec(b, key, n, who) -> Vector{Float64}

Lee un vector numérico obligatorio de `n` componentes.

# Argumentos
- `b::AbstractDict`: bloque leído del JSON.
- `key::AbstractString`: nombre del parámetro.
- `n::Int`: número de componentes esperado.
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.

# Devuelve
- `Vector{Float64}` de longitud `n`.

# Excepciones
- `ErrorException`: si falta, no es un vector de `n` números o tiene valores no finitos.
"""
function getvec(b::AbstractDict, key::AbstractString, n::Int, who::AbstractString)
    haskey(b, key) || error("$who: missing '$key'.")
    v = b[key]
    (v isa AbstractVector && length(v) == n && all(x -> x isa Real && !(x isa Bool), v)) ||
        error("$who: '$key' must be a list of $n numbers (got $(repr(v))).")
    out = Float64.(v)
    all(isfinite, out) || error("$who: '$key' must have finite values.")
    return out
end

"""
    getunit(b, key, who) -> SVector{3,Float64}

Lee un vector de 3 componentes y lo normaliza.

# Argumentos
- `b::AbstractDict`: bloque leído del JSON.
- `key::AbstractString`: nombre del parámetro.
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.

# Devuelve
- `SVector{3,Float64}` unitario.

# Excepciones
- `ErrorException`: si falta, no tiene 3 componentes o su norma es cero.
"""
function getunit(b::AbstractDict, key::AbstractString, who::AbstractString)
    v = Vec3(getvec(b, key, 3, who))
    n = norm(v)
    n > 0 || error("$who: '$key' has zero norm.")
    return v / n
end

# -----------------------------------------------------------------------------
# Lectura del bloque "attitude"
# -----------------------------------------------------------------------------

"""
    parse_attitude(block, name) -> AttitudeSpec

Lee y valida un bloque `attitude` con el modo indicado en su campo `type`.

# Argumentos
- `block::AbstractDict`: bloque `attitude` (de JSON3 o de un `Dict`).
- `name::AbstractString`: nombre del objetivo, para los mensajes de error.

# Devuelve
- Un subtipo de `AttitudeSpec` (`SpinSpec`, `RandomSpec`, `FixedSpec`, `PointingSpec`...).

# Excepciones
- `ErrorException`: si el bloque no es un diccionario, no tiene `type` o el modo no existe.
- `ErrorException`: cualquier error de validación del modo elegido.

# Ejemplo
```julia
parse_attitude(Dict("type" => "spin", "w_deg_s" => 2, "w_dir" => [0, 0, 1]), "deb")
# SpinSpec(...)
```
"""
function parse_attitude(block, name::AbstractString)
    b = todict(block)
    b isa Dict{String,Any} || error("Objective '$name': 'attitude' must be a block {...} (got $(repr(block))).")
    haskey(b, "type") || error("Objective '$name': 'attitude' block without 'type'.")

    # El modo se busca en el registro; si no existe se listan los disponibles
    mode = String(b["type"])
    haskey(ATTITUDES, mode) || error(
        "Objective '$name': unknown attitude type '$mode'. Available: " * join(sort!(collect(keys(ATTITUDES))), ", "))

    return ATTITUDES[mode](b, String(name))
end

"""
    legacy_to_spin(w_deg_s, w_dir=nothing) -> Dict{String,Any}

Traduce las claves antiguas `w_deg_s` / `w_dir` a un bloque `attitude` de tipo `spin`.

Por qué existe: `attitude_block` la usa para mostrar, en el mensaje de error, el
bloque que sustituye a las claves antiguas de un JSON.

# Argumentos
- `w_deg_s::Real`: velocidad angular (°/s).
- `w_dir`: eje de giro en ejes cuerpo, o `nothing` si no se dio.

# Devuelve
- `Dict{String,Any}` con `type = "spin"`, `w_deg_s` y, si se dio, `w_dir`.

# Excepciones
- Ninguna (la validación la hace `parse_attitude`).

# Ejemplo
```julia
legacy_to_spin(0)                  # Dict("type" => "spin", "w_deg_s" => 0)
legacy_to_spin(30, [1, 0, 0])      # Dict("type" => "spin", "w_deg_s" => 30, "w_dir" => [1, 0, 0])
```
"""
function legacy_to_spin(w_deg_s, w_dir = nothing)
    block = Dict{String,Any}("type" => "spin", "w_deg_s" => w_deg_s)
    w_dir === nothing || (block["w_dir"] = collect(w_dir))
    return block
end

"""
    attitude_block(entry, name) -> Dict{String,Any}

Devuelve el bloque `attitude` de un objetivo (del JSON del simulador o de un
grupo de la spec) tras comprobar que está bien definido. Las claves antiguas
`w_deg_s` / `w_dir` ya no se aceptan: dan un error que muestra el bloque
equivalente para copiarlo.

Por qué existe: el simulador y el generador aplican exactamente las mismas reglas
y los mismos mensajes de error.

# Argumentos
- `entry::AbstractDict`: entrada del objetivo (o grupo de objetivos de la spec).
- `name::AbstractString`: nombre del objetivo, para los mensajes de error.

# Devuelve
- `Dict{String,Any}`: copia del bloque `attitude` (sin validar sus parámetros).

# Excepciones
- `ErrorException`: si falta el bloque `attitude`, si no es un bloque `{...}` o si
  aparecen `w_deg_s` / `w_dir` (el mensaje indica cómo convertirlas).

# Ejemplo
```julia
attitude_block(Dict("attitude" => Dict("type" => "spin", "w_deg_s" => 0)), "deb")
# Dict("type" => "spin", "w_deg_s" => 0)
attitude_block(Dict("w_deg_s" => 0, "w_dir" => [1, 0, 0]), "deb")
# error: ... Replace them with:  "attitude": {"type": "spin", "w_deg_s": 0, "w_dir": [1, 0, 0]}
```
"""
function attitude_block(entry, name::AbstractString)
    e = todict(entry)

    # Formato antiguo: se indica el bloque exacto que hay que escribir
    if haskey(e, "w_deg_s") || haskey(e, "w_dir")
        haskey(e, "w_deg_s") || error("Objective '$name': 'w_dir' without 'w_deg_s'; use an 'attitude' block.")
        block = legacy_to_spin(e["w_deg_s"], get(e, "w_dir", nothing))
        # Valores escritos como en JSON: textos entre comillas y vectores como [1, 0, 0]
        fmt(x) = x isa AbstractVector ? "[" * join(fmt.(x), ", ") * "]" :
                 x isa AbstractString ? "\"$x\"" : string(x)
        txt = "{" * join(["\"$k\": $(fmt(block[k]))" for k in ("type", "w_deg_s", "w_dir") if haskey(block, k)], ", ") * "}"
        error("Objective '$name': 'w_deg_s'/'w_dir' are no longer accepted. " *
              "Replace them with:  \"attitude\": $txt")
    end

    haskey(e, "attitude") || error(
        "Objective '$name': missing 'attitude' block (e.g. \"attitude\": {\"type\": \"spin\", \"w_deg_s\": 0}).")
    e["attitude"] isa Dict{String,Any} || error("Objective '$name': 'attitude' must be a block {...}.")
    return e["attitude"]
end

"""
    parse_member_attitude(entry, name) -> AttitudeSpec

Lee la actitud de un objetivo del JSON del simulador: comprueba el bloque con
`attitude_block` y lo valida con `parse_attitude`.

# Argumentos
- `entry::AbstractDict`: entrada del objetivo en la sección `objectives`.
- `name::AbstractString`: nombre del objetivo, para los mensajes de error.

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
function parse_member_attitude(entry, name::AbstractString)
    return parse_attitude(attitude_block(entry, name), name)
end

# -----------------------------------------------------------------------------
# Cálculo de la actitud
# -----------------------------------------------------------------------------

"""
    compute_attitude(spec, traj, sun_dir) -> AttitudeData

Calcula la historia de la DCM cuerpo -> ECI de un objetivo. Cada modo añade su
propio método; este es el de reserva, que solo da error.

# Argumentos
- `spec::AttitudeSpec`: actitud leída con `parse_attitude`.
- `traj::TrajectoryData`: trayectoria del objetivo (r, v en ECI, km y km/s; t en s).
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
