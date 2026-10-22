# =============================================================================
# cw_state.jl
#
# Plugin de disposición: estado relativo inicial dado a mano, en el sistema LVLH
# del chief. Registra la disposición "cw_state". Es el caso general de CW: cw_ellipse
# y pco generan estos mismos estados a partir de una geometría.
#
#   "layout": {"type": "cw_state", "states": [
#       {"x_km": 0.5, "y_km": 0.0, "z_km": 0.0, "vx_km_s": 0.0, "vy_km_s": -0.00108, "vz_km_s": 0.0},
#       {"x_km": 0.0, "y_km": 2.0, "z_km": 0.3, "vx_km_s": 0.0, "vy_km_s": 0.0, "vz_km_s": 0.0}
#   ]}
#
# Ejes LVLH: x radial (afuera), y a lo largo de la traza, z normal al plano. Las seis
# claves de cada estado son obligatorias (una velocidad omitida no es "cero" a
# efectos de deriva: vy distinto de -2·n·x produce deriva secular).
# Opcional: "include_chief": true añade el chief como primer miembro.
#
# Funciones: cw_state_layout
# =============================================================================

# Claves obligatorias de cada estado relativo
const CW_STATE_KEYS = ("x_km", "y_km", "z_km", "vx_km_s", "vy_km_s", "vz_km_s")

"""
    cw_state_layout(chief, layout, ctx) -> Vector{Dict{String,Any}}

Genera un miembro por cada estado relativo LVLH de la lista `states`.

# Argumentos
- `chief::AbstractDict`: elementos keplerianos del chief (de `build_orbit`).
- `layout::AbstractDict`: bloque `layout` con la lista `states` (y `include_chief` opcional).
- `ctx`: contexto del Builder.

# Devuelve
- `Vector` de diccionarios de elementos, uno por estado y en el mismo orden.

# Excepciones
- `ErrorException`: si falta `states`, está vacía, o un estado no tiene exactamente las seis claves.
- `ErrorException`: si un estado no produce una órbita ligada.

# Ejemplo
```julia
layout = Dict("type" => "cw_state", "states" => [
    Dict("x_km" => 0.0, "y_km" => 1.0, "z_km" => 0.0,
         "vx_km_s" => 0.0, "vy_km_s" => 0.0, "vz_km_s" => 0.0)])
cw_state_layout(chief, layout, ctx)      # un miembro 1 km por delante del chief
```
"""
function cw_state_layout(chief::AbstractDict, layout::AbstractDict, ctx)
    check_keys(layout, ("type", "states", "include_chief"); who = "cw_state")

    haskey(layout, "states") || error("'cw_state': missing required parameter 'states'.")
    raw = layout["states"]
    isempty(raw) && error("'cw_state': 'states' is empty.")

    states = Vector{Float64}[]
    for (k, s) in enumerate(raw)
        check_keys(s, CW_STATE_KEYS; who = "cw_state (state $k)")
        push!(states, [getreq(s, key) for key in CW_STATE_KEYS])    # orden [x y z vx vy vz]
    end

    return cw_members(chief, states, layout, ctx; who = "cw_state")
end

register_layout!("cw_state", cw_state_layout)
