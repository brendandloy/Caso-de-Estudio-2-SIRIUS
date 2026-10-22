# =============================================================================
# cw_ellipse.jl
#
# Plugin de disposición: N satélites sobre una misma elipse relativa de CW
# alrededor del chief, repartidos uniformemente en fase. Registra "cw_ellipse".
#
#   "layout": {"type": "cw_ellipse", "n": 3, "x_amplitude_km": 0.5,
#              "y_center_km": 0.0, "z_amplitude_km": 0.3, "phase_deg": 0}
#
# Con θ = fase del miembro (θ_k = phase_deg + 360·k/n; en el tiempo θ = n_medio·t + θ_k):
#   x =  A·cos θ                  radial
#   y =  y_c − 2A·sin θ           a lo largo de la traza (semieje doble que el radial)
#   z =  B·cos(θ − z_phase)       fuera de plano
# con A = x_amplitude_km, B = z_amplitude_km, y_c = y_center_km. La velocidad inicial
# cumple ẏ0 = −2·n·x0 (misma energía: sin deriva secular entre los miembros y el chief).
#
# Obligatorios: n, x_amplitude_km, y_center_km, z_amplitude_km.
# Opcionales: phase_deg (0), z_phase_deg (0), include_chief (false).
#
# Funciones: cw_ellipse_layout
# =============================================================================

"""
    cw_ellipse_layout(chief, layout, ctx) -> Vector{Dict{String,Any}}

Genera `n` miembros sobre la elipse relativa CW definida por las amplitudes y el
centro, y los convierte a elementos keplerianos.

# Argumentos
- `chief::AbstractDict`: elementos keplerianos del chief (casi circular).
- `layout::AbstractDict`: bloque `layout` con `n`, `x_amplitude_km`, `y_center_km`,
  `z_amplitude_km` y, opcionalmente, `phase_deg`, `z_phase_deg`, `include_chief`.
- `ctx`: contexto del Builder (usa `ctx.mu`).

# Devuelve
- `Vector` de `n` diccionarios de elementos (uno más si `include_chief`).

# Excepciones
- `ErrorException`: si aparece una clave no admitida o falta una obligatoria.
- `ErrorException`: si `n` no es un entero >= 1, `x_amplitude_km <= 0` o `z_amplitude_km < 0`.

# Ejemplo
```julia
layout = Dict("type" => "cw_ellipse", "n" => 2, "x_amplitude_km" => 0.5,
              "y_center_km" => 0.0, "z_amplitude_km" => 0.3)
cw_ellipse_layout(chief, layout, ctx)      # 2 miembros en fases 0 y 180 grados de la elipse
```
"""
function cw_ellipse_layout(chief::AbstractDict, layout::AbstractDict, ctx)
    check_keys(layout, ("type", "n", "x_amplitude_km", "y_center_km", "z_amplitude_km",
                        "phase_deg", "z_phase_deg", "include_chief"); who = "cw_ellipse")

    nf = getreq(layout, "n")
    (isinteger(nf) && nf >= 1) || error("'cw_ellipse': 'n' must be an integer >= 1 (got $nf).")
    n = Int(nf)

    A  = getreq(layout, "x_amplitude_km")
    yc = getreq(layout, "y_center_km")
    B  = getreq(layout, "z_amplitude_km")
    A > 0  || error("'cw_ellipse': 'x_amplitude_km' must be > 0 (got $A).")
    B >= 0 || error("'cw_ellipse': 'z_amplitude_km' must be >= 0 (got $B).")

    phase = getopt(layout, "phase_deg", 0.0)
    zph   = deg2rad(getopt(layout, "z_phase_deg", 0.0))
    nm    = sqrt(ctx.mu / chief["SMA"]^3)                       # movimiento medio del chief (rad/s)

    states = Vector{Float64}[]
    for k in 0:n-1
        θ = deg2rad(phase + 360.0 * k / n)                      # fase del miembro k en la elipse
        push!(states, [
            A * cos(θ),                                         # x
            yc - 2 * A * sin(θ),                                # y
            B * cos(θ - zph),                                   # z
            -A * nm * sin(θ),                                   # ẋ
            -2 * A * nm * cos(θ),                               # ẏ = -2·n·x (sin deriva)
            -B * nm * sin(θ - zph)                              # ż
        ])
    end

    return cw_members(chief, states, layout, ctx; who = "cw_ellipse")
end

register_layout!("cw_ellipse", cw_ellipse_layout)
