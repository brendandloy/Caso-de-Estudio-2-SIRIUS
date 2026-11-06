# =============================================================================
# JSON_utils.jl
#
# Módulo JSONUtils: lectores y validadores comunes de los bloques del JSON (spec y
# JSON completo del simulador). Los usan Attitude, Geometry y Scenario, para que
# todos apliquen las mismas reglas y den los mismos mensajes de error.
#
# Funciones:
#   todict           : objeto JSON3 -> Dict{String,Any} / Vector{Any}, recursivo
#   check_keys       : rechaza claves no admitidas (y las calculadas por el tipo)
#   getnum           : número obligatorio y finito
#   getvec           : vector numérico obligatorio de n componentes
#   getunit          : vector de 3 componentes normalizado
#   get_or_draw_seed : semilla del bloque o una nueva sorteada (con aviso)
#
# Dependencias: ninguna del proyecto; se incluye antes que Attitude.
# =============================================================================

module JSONUtils

export todict, check_keys, getnum, getvec, getunit, get_or_draw_seed

using LinearAlgebra
using Random
using StaticArrays

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
    check_keys(b, allowed; who, computed=()) -> Nothing

Comprueba que un bloque del JSON no tiene claves fuera de las admitidas.

Por qué existe: una errata (p. ej. `w_deg` en vez de `w_deg_s`, o `"RAAm"`) se
ignoraría en silencio y la simulación usaría otro valor sin avisar. Con `computed`
impide además fijar a mano un valor que el tipo calcula solo (p. ej. `INC` en una
órbita SSO). La usan Attitude, Geometry y Scenario.

# Argumentos
- `b::AbstractDict`: bloque leído del JSON (claves `String` o `Symbol`).
- `allowed`: colección de claves admitidas.
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.
- `computed`: claves que el tipo de órbita calcula por sí mismo; si aparecen, el error lo explica.

# Devuelve
- `nothing` si todas las claves son válidas.

# Excepciones
- `ErrorException`: si aparece una clave de `computed` (no se puede fijar a mano).
- `ErrorException`: si aparece cualquier otra clave que no esté en `allowed` (las lista todas).

# Ejemplo
```julia
check_keys(Dict("type" => "spin", "w_deg" => 1), ("type", "w_deg_s"); who = "deb")
# error: deb: unknown key(s) w_deg. Allowed: type, w_deg_s
```
"""
function check_keys(b::AbstractDict, allowed; who::AbstractString, computed = ())
    allowed_set = Set(String.(collect(allowed)))
    bad = sort!([String(k) for k in keys(b) if !(String(k) in allowed_set)])   # claves no admitidas
    isempty(bad) && return nothing

    # Claves que el tipo calcula solo: error específico, con la alternativa
    fixed = filter(in(Set(String.(collect(computed)))), bad)
    isempty(fixed) || error("$who: $(join(fixed, ", ")) fixed/computed by this orbit type and cannot be set. " *
                            "Use type 'keplerian' if you need to define it manually.")

    error("$who: unknown key(s) $(join(bad, ", ")). Allowed: $(join(sort!(collect(allowed_set)), ", "))")
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
    v = SVector{3,Float64}(getvec(b, key, 3, who))
    n = norm(v)
    n > 0 || error("$who: '$key' has zero norm.")
    return v / n
end

"""
    get_or_draw_seed(b, who; hint) -> Int

Lee la semilla `"seed"` de un bloque o, si no la tiene, sortea una nueva y la
imprime para que se pueda repetir el resultado.

Por qué existe: la actitud `random`, la forma de cada objetivo y la spec completa
siguen la misma regla con la semilla; así se valida y se avisa igual en todos.

# Argumentos
- `b::AbstractDict`: bloque leído del JSON.
- `who::AbstractString`: texto que identifica el bloque en los mensajes.
- `hint::AbstractString`: dónde añadir la semilla, para el aviso (p. ej. `"to its 'shape' block"`).

# Devuelve
- `Int`: la semilla leída o sorteada (≥ 0).

# Excepciones
- `ErrorException`: si `"seed"` existe pero no es un entero ≥ 0.

# Ejemplo
```julia
get_or_draw_seed(Dict("seed" => 42), "deb"; hint = "to its 'shape' block")   # 42
```
"""
function get_or_draw_seed(b::AbstractDict, who::AbstractString; hint::AbstractString)
    if haskey(b, "seed")
        s = b["seed"]
        (s isa Real && !(s isa Bool) && isinteger(s) && s >= 0) ||
            error("$who: 'seed' must be an integer >= 0 (got $(repr(s))).")
        return Int(s)
    end

    # Sin semilla: una nueva de la entropía del sistema, impresa para poder repetir
    seed = Int(rand(RandomDevice(), UInt32))
    println("$who: no 'seed'; using seed = $seed (add \"seed\": $seed $hint to repeat it).")
    return seed
end

end # module JSONUtils
