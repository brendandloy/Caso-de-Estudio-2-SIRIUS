# =============================================================================
# pco.jl
#
# Plugin de disposición: órbita circular proyectada (PCO). Registra "pco".
# Es cw_ellipse con la amplitud en z elegida para que la proyección sobre el plano
# y–z (visto desde el eje radial) sea una circunferencia.
#
#   "layout": {"type": "pco", "n": 3, "radius_km": 1.0, "y_center_km": 0.0}
#
# radius_km es el radio de la circunferencia proyectada: A = radius/2 (radial) y
# B = radius (z). La distancia 3D no es constante: oscila entre `radius` y 1.118·radius.
#
# Obligatorios: n, radius_km, y_center_km. Opcionales: phase_deg (0), include_chief (false).
#
# Funciones: pco_layout
# =============================================================================

"""
    pco_layout(chief, layout, ctx) -> Vector{Dict{String,Any}}

Genera `n` miembros repartidos uniformemente sobre una PCO de radio `radius_km`,
delegando en `cw_ellipse_layout` con `x_amplitude_km = radius/2` y
`z_amplitude_km = radius`.

# Argumentos
- `chief::AbstractDict`: elementos keplerianos del chief (casi circular).
- `layout::AbstractDict`: bloque `layout` con `n`, `radius_km`, `y_center_km` y,
  opcionalmente, `phase_deg`, `include_chief`.
- `ctx`: contexto del Builder.

# Devuelve
- `Vector` de `n` diccionarios de elementos (uno más si `include_chief`).

# Excepciones
- `ErrorException`: si aparece una clave no admitida o falta una obligatoria.
- `ErrorException`: si `radius_km <= 0` o cualquier error de `cw_ellipse_layout`.

# Ejemplo
```julia
layout = Dict("type" => "pco", "n" => 3, "radius_km" => 1.0, "y_center_km" => 0.0)
pco_layout(chief, layout, ctx)      # 3 miembros a 120 grados sobre una circunferencia de 1 km
```
"""
function pco_layout(chief::AbstractDict, layout::AbstractDict, ctx)
    check_keys(layout, ("type", "n", "radius_km", "y_center_km", "phase_deg", "include_chief");
               who = "pco")

    radius = getreq(layout, "radius_km")
    radius > 0 || error("'pco': 'radius_km' must be > 0 (got $radius).")

    # Se traduce a la elipse equivalente y se reutiliza cw_ellipse
    ellipse = Dict{String,Any}(
        "type"           => "cw_ellipse",
        "n"              => getreq(layout, "n"),
        "x_amplitude_km" => radius / 2,
        "y_center_km"    => getreq(layout, "y_center_km"),
        "z_amplitude_km" => radius,
        "phase_deg"      => getopt(layout, "phase_deg", 0.0),
        "include_chief"  => get(layout, "include_chief", false)
    )
    return cw_ellipse_layout(chief, ellipse, ctx)
end

register_layout!("pco", pco_layout)
