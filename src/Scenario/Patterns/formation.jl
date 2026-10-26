# =============================================================================
# formation.jl
#
# Patrón de formación: varios satélites definidos respecto a una órbita jefe.
# Registra el patrón "formation".
#
#   "cam": {
#     "pattern":  "formation",
#     "chief":    {"type": "LEO", "altitude_km": 620, "INC": 90, "RAAN": 0, "TA": 0},
#     "layout":   {"type": "along_track", "n": 3, "spacing_km": 10},
#     "pointing": {"type": "fixed", "dir": [0, 1, 0]}
#   }
#
# - "chief": bloque "orbit" de cualquier tipo registrado (LEO, SSO, GEO, keplerian...).
# - "layout": cómo se reparten los miembros (ver carpeta Layouts/).
# Los miembros se llaman cam1..camN (los nombra el Builder). El resto de campos del
# grupo (pointing, sigmas, attitude...) se aplican igual a todos los miembros.
#
# Funciones: formation_pattern
# =============================================================================

"""
    formation_pattern(group, ctx) -> Vector{Dict{String,Any}}

Construye la órbita jefe y delega en la disposición indicada para generar los
miembros de la formación.

Por qué existe: separa "qué órbita" (cualquier plugin de Orbits/) de "cómo se
reparten los satélites" (cualquier plugin de Layouts/); se combinan libremente.

# Argumentos
- `group::AbstractDict`: grupo de la spec con los bloques `chief` y `layout`.
- `ctx`: contexto del Builder.

# Devuelve
- `Vector` de diccionarios de elementos keplerianos, uno por miembro.

# Excepciones
- `ErrorException`: si falta el bloque `chief` o `layout`, o `layout` no tiene `type`.
- `ErrorException`: si el tipo de disposición no está registrado.
- `ErrorException`: cualquier error del plugin de órbita o de la disposición.

# Ejemplo
```julia
group = Dict(
    "chief"  => Dict("type" => "GEO", "longitude_deg" => 0.0),
    "layout" => Dict("type" => "along_track", "n" => 2, "spacing_km" => 50.0))
formation_pattern(group, ctx)          # 2 miembros GEO separados 50 km en la traza
```
"""
function formation_pattern(group::AbstractDict, ctx)
    haskey(group, "chief")  || error("Pattern 'formation' needs a 'chief' orbit block.")
    haskey(group, "layout") || error("Pattern 'formation' needs a 'layout' block.")

    layout = group["layout"]
    haskey(layout, "type") || error("'layout' block without 'type'.")

    chief = build_orbit(group["chief"], ctx)
    return lookup(LAYOUTS, "layout", String(layout["type"]))(chief, layout, ctx)
end

register_pattern!("formation", formation_pattern)
