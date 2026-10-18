# =============================================================================
# single.jl
#
# Patrón más simple: un único objeto con una órbita. Registra el patrón "single".
#
#   "cam1": {"pattern": "single", "orbit": {"type": "leo", "altitude_km": 621.86, "INC": 90,
#                                           "RAAN": 0, "TA": 0}, ...}
#
# Si se omite "pattern", el Builder usa "single" por defecto.
#
# Funciones: single_pattern
# =============================================================================

"""
    single_pattern(group, ctx) -> Vector{Dict{String,Any}}

Devuelve un único miembro cuya órbita se construye con el bloque `orbit` del grupo.

Por qué existe: el Builder siempre espera un patrón que devuelva una lista de
miembros. Este es el caso trivial de una lista con un solo elemento, y sirve de
modelo para patrones que devuelven varios (formación, constelación).

# Argumentos
- `group::AbstractDict`: grupo de la spec; debe contener un bloque `orbit`.
- `ctx`: contexto del Builder.

# Devuelve
- `Vector` con un solo diccionario de elementos keplerianos.

# Excepciones
- `ErrorException`: si el grupo no tiene bloque `orbit`.
- `ErrorException`: cualquier error del plugin de órbita elegido.

# Ejemplo
```julia
group = Dict("orbit" => Dict("type" => "geo", "longitude_deg" => 0.0))
single_pattern(group, ctx)        # [Dict(...)] con un solo miembro
```
"""
function single_pattern(group::AbstractDict, ctx)
    haskey(group, "orbit") || error("Pattern 'single' needs an 'orbit' block.")
    return [build_orbit(group["orbit"], ctx)]
end

register_pattern!("single", single_pattern)
