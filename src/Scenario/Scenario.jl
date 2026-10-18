module Scenario

# =============================================================================
# Módulo Scenario
#
# Generador de JSON de escenario para SIRIUS. Convierte una spec corta
# (config/specs/*.json) en un JSON completo, con el mismo esquema que ya lee
# Environment.load_config.
#
# Scripts que contiene:
#   Scenario.jl        : este archivo; registros, helpers comunes y autocarga de plugins
#   Builder.jl         : valores por defecto, validación, construcción y escritura del JSON
#   Orbital_utils.jl   : helpers de física orbital compartidos por los plugins de órbita
#   Orbits/            : keplerian.jl, leo.jl, sso.jl, geo.jl
#   Patterns/          : single.jl
#   Pointing/          : fixed.jl
#
# Hay tres tipos de "plugin", cada uno en su carpeta y con un contrato fijo:
#   Orbits/    f(params, ctx)          -> Dict con SMA, ECC, INC, RAAN, AOP, TA
#   Patterns/  f(group,  ctx)          -> Vector{Dict} de elementos (miembros)
#   Pointing/  f(spec, member, ctx)    -> Vector{Float64}(3) = obs_dir
#
# Para añadir un tipo nuevo basta con crear un archivo en la carpeta que
# corresponda y terminarlo con register_orbit! / register_pattern! /
# register_pointing!. No hay que tocar nada más: este módulo incluye
# automáticamente todos los .jl de esas carpetas.
#
# Dependencias: Constants (μ, R, ω_earth) y Transformations (theta0_from_utc),
# que deben estar cargados antes (ver scripts/generate_config.jl).
# =============================================================================

export generate_config, generate_config_file, write_config
export register_orbit!, register_pattern!, register_pointing!

using ..Constants          # μ, R (km), ω_earth (rad/s)
using ..Transformations: theta0_from_utc   # misma rotación terrestre (ERA) que el simulador
using Dates
using JSON3
using LinearAlgebra

# -----------------------------------------------------------------------------
# Registros
# -----------------------------------------------------------------------------

# Un registro por tipo de plugin: nombre del tipo -> función que lo implementa
const ORBITS    = Dict{String,Function}()
const PATTERNS  = Dict{String,Function}()
const POINTINGS = Dict{String,Function}()

# Registra una función de órbita bajo un nombre; devuelve nothing
register_orbit!(name::AbstractString, f::Function)    = (ORBITS[String(name)]    = f; nothing)

# Registra una función de patrón bajo un nombre; devuelve nothing
register_pattern!(name::AbstractString, f::Function)  = (PATTERNS[String(name)]  = f; nothing)

# Registra una función de apuntado bajo un nombre; devuelve nothing
register_pointing!(name::AbstractString, f::Function) = (POINTINGS[String(name)] = f; nothing)

"""
    lookup(registry, kind, name) -> Function

Busca en un registro la función asociada a un nombre de tipo.

Por qué existe: es lo que permite añadir tipos nuevos sin editar el código
existente. El Builder no conoce ningún tipo concreto; solo pregunta al
registro por el nombre que aparece en la spec.

# Argumentos
- `registry::AbstractDict`: registro donde buscar (ORBITS, PATTERNS o POINTINGS).
- `kind::AbstractString`: texto que describe el registro, solo para el mensaje de error.
- `name::AbstractString`: nombre del tipo pedido en la spec.

# Devuelve
- `Function`: la función registrada con ese nombre.

# Excepciones
- `ErrorException`: si el nombre no está registrado; el mensaje lista los disponibles.

# Ejemplo
```julia
f = lookup(ORBITS, "orbit", "GEO")        # devuelve geo_elements
lookup(ORBITS, "orbit", "geo")            # error: Unknown orbit type 'geo'. Available: GEO, LEO, SSO, keplerian
```
"""
function lookup(registry::AbstractDict, kind::AbstractString, name::AbstractString)
    # Si no existe, se informa de los tipos que sí hay para detectar erratas
    haskey(registry, name) || error(
        "Unknown $kind type '$name'. Available: " * join(sort!(collect(keys(registry))), ", ")
    )
    return registry[name]
end

# -----------------------------------------------------------------------------
# Helpers para los plugins
# -----------------------------------------------------------------------------

# Lee un parámetro obligatorio como Float64; error claro si falta
getreq(p::AbstractDict, key::AbstractString) =
    haskey(p, key) ? Float64(p[key]) : error("Missing required parameter '$key'.")

# Lee un parámetro opcional como Float64; devuelve `default` si falta
getopt(p::AbstractDict, key::AbstractString, default::Real) =
    haskey(p, key) ? Float64(p[key]) : Float64(default)

"""
    build_orbit(spec, ctx) -> Dict{String,Any}

Construye los seis elementos keplerianos de un bloque `orbit`, delegando en
el plugin de órbita registrado con su campo `type`.

Por qué existe: los patrones (single, formation, walker...) la usan para
obtener la órbita base sin saber qué tipos de órbita hay. Así una formación
funciona con cualquier tipo, incluidos los que añadas en el futuro.

# Argumentos
- `spec::AbstractDict`: bloque `orbit` de la spec (debe tener `type`).
- `ctx`: contexto del Builder (época, constantes, etc.).

# Devuelve
- `Dict{String,Any}` con `SMA, ECC, INC, RAAN, AOP, TA`.

# Excepciones
- `ErrorException`: si el bloque no tiene `type`.
- `ErrorException`: si el tipo no está registrado o el plugin rechaza los parámetros.

# Ejemplo
```julia
orbit = Dict("type" => "GEO", "longitude_deg" => 0.0)
el = build_orbit(orbit, ctx)              # Dict con SMA ≈ 42164.17, INC = 0, TA ≈ 280.33
```
"""
function build_orbit(spec::AbstractDict, ctx)
    haskey(spec, "type") || error("Orbit block without 'type'.")
    return lookup(ORBITS, "orbit", String(spec["type"]))(spec, ctx)
end

# -----------------------------------------------------------------------------
# Núcleo + autocarga de plugins
# -----------------------------------------------------------------------------

include(joinpath(@__DIR__, "Builder.jl"))
include(joinpath(@__DIR__, "Orbital_utils.jl"))   # helpers compartidos por los plugins de órbita

# Se incluyen todos los .jl de cada carpeta de plugins, en orden alfabético
for sub in ("Orbits", "Patterns", "Pointing")
    dir = joinpath(@__DIR__, sub)
    isdir(dir) || continue

    # Si el módulo se precompila, esto permite detectar archivos nuevos; si falla, no importa
    try
        include_dependency(dir)
    catch
    end

    for f in sort!(filter(endswith(".jl"), readdir(dir)))
        include(joinpath(dir, f))
    end
end

end # module Scenario