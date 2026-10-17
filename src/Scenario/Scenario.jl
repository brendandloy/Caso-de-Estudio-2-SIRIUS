module Scenario

# =============================================================================
# Generador de JSON de escenario para SIRIUS
#
# spec corto (config/specs/*.json)  ->  JSON completo (mismo esquema que ya
# lee Environment.load_config).
#
# Tres tipos de "plugin", cada uno en su carpeta y con un contrato fijo:
#   Orbits/    f(params, ctx)          -> Dict con SMA, ECC, INC, RAAN, AOP, TA
#   Patterns/  f(group,  ctx)          -> Vector{Dict} de elementos (miembros)
#   Pointing/  f(spec, member, ctx)    -> Vector{Float64}(3) = obs_dir
#
# Para añadir un tipo nuevo: crear un archivo en la carpeta correspondiente
# y terminarlo con register_orbit! / register_pattern! / register_pointing!.
# No hay que tocar nada más.
# =============================================================================

export generate_config, generate_config_file, write_config
export register_orbit!, register_pattern!, register_pointing!

using ..Constants          # μ, R (km)
using Dates
using JSON3
using LinearAlgebra

# -----------------------------------------------------------------------------
# Registros
# -----------------------------------------------------------------------------

# Diccionarios globales que actúan como registros para almacenar funciones asociadas a sus nombres (claves)
const ORBITS    = Dict{String,Function}()
const PATTERNS  = Dict{String,Function}()
const POINTINGS = Dict{String,Function}()

# Registra una función de órbita en el diccionario global ORBITS y devuelve nothing
register_orbit!(name::AbstractString, f::Function)    = (ORBITS[String(name)]    = f; nothing)

# Registra una función de patrón en el diccionario global PATTERNS y devuelve nothing
register_pattern!(name::AbstractString, f::Function)  = (PATTERNS[String(name)]  = f; nothing)

# Registra una función de apuntamiento en el diccionario global POINTINGS y devuelve nothing
register_pointing!(name::AbstractString, f::Function) = (POINTINGS[String(name)] = f; nothing)

function lookup(registry::AbstractDict, kind::AbstractString, name::AbstractString)
    """
    Busca y devuelve una función registrada en un diccionario dado por su nombre.

    # Argumentos
    - `registry::AbstractDict`: Diccionario que actúa como registro (ej. `ORBITS`, `PATTERNS`, `POINTINGS`).
    - `kind::AbstractString`: Tipo o categoría del registro, utilizado para personalizar el mensaje de error.
    - `name::AbstractString`: Nombre o clave de la función que se desea obtener.

    # Devuelve
    - `Function`: La función asociada a la clave `name` dentro del registro.

    # Excepciones
    - `ErrorException`: Se lanza si el nombre especificado no existe en el registro, indicando las claves disponibles.
    """

    # Comprueba si la clave existe; si no, lanza un error listando las claves disponibles ordenadas alfabéticamente
    haskey(registry, name) || error(
        "Unknown $kind type '$name'. Available: " * join(sort!(collect(keys(registry))), ", ")
    )

    # Devuelve la función correspondiente a la clave solicitada
    return registry[name]
end

# -----------------------------------------------------------------------------
# Helpers para los plugins
# -----------------------------------------------------------------------------

# Extraer un parámetro requerido de un diccionario, o lanzar un error si no existe.
getreq(p::AbstractDict, key::AbstractString) =
    haskey(p, key) ? Float64(p[key]) : error("Missing required parameter '$key'.")

# Extraer un parámetro opcional de un diccionario, o devolver un valor por defecto si no existe.
getopt(p::AbstractDict, key::AbstractString, default::Real) =
    haskey(p, key) ? Float64(p[key]) : Float64(default)

function build_orbit(spec::AbstractDict, ctx)
    """
    Instancia dinámicamente una órbita según el tipo especificado en su configuración.

    # Argumentos
    - `spec::AbstractDict`: Diccionario de configuración de la órbita (debe incluir la clave `"type"`).
    - `ctx`: Contexto de simulación o ejecución necesario para la construcción de la órbita.

    # Devuelve
    - Instancia de la órbita generada por la función constructora correspondiente.

    # Excepciones
    - `ErrorException`: Se lanza si `spec` no contiene la clave `"type"`.
    - `ErrorException`: Se lanza (vía `lookup`) si el tipo especificado no está registrado en `ORBITS`.
    """

    # Verifica que el diccionario contenga la clave "type"; si no existe, lanza un error
    haskey(spec, "type") || error("Orbit block without 'type'.")

    # Recupera la función constructora desde el registro ORBITS a partir del tipo indicado
    funcion_constructora = lookup(ORBITS, "orbit", String(spec["type"]))

    # Ejecuta el constructor pasándole el diccionario de especificaciones y el contexto
    resultado = funcion_constructora(spec, ctx)

    # Devuelve el objeto u órbita instanciada
    return resultado
end

# -----------------------------------------------------------------------------
# Núcleo + autocarga de plugins
# -----------------------------------------------------------------------------

# Carga el código base/núcleo ubicado en la misma carpeta
include(joinpath(@__DIR__, "Builder.jl"))

# Itera sobre las subcarpetas de plugins conocidas
for sub in ("Orbits", "Patterns", "Pointing")
    dir = joinpath(@__DIR__, sub)
    isdir(dir) || continue  # Salta la carpeta si no existe

    # Registra la carpeta como dependencia para detectar cambios durante la precompilación
    try
        include_dependency(dir)
    catch
    end

    # Lee, ordena alfabéticamente e incluye todos los archivos .jl de la subcarpeta
    for f in sort!(filter(endswith(".jl"), readdir(dir)))
        include(joinpath(dir, f))
    end
end

end # module Scenario