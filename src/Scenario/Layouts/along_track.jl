# =============================================================================
# along_track.jl
#
# Plugin de disposición: miembros en la MISMA órbita, separados a lo largo de
# la traza. Registra la disposición "along_track" (para el patrón "formation").
#
#   "layout": {"type": "along_track", "n": 3, "spacing_km": 10}
#
# El miembro k (k = 0..n-1) va k·spacing_km por delante del jefe (el miembro 1
# es el jefe). spacing_km < 0 los coloca detrás.
#
# Física: misma SMA, ECC, INC, RAAN y AOP; solo cambia la anomalía media,
# ΔM = spacing / SMA. Con la misma SMA todos tienen el mismo periodo, así que la
# separación se mantiene en el modelo kepleriano (con J2 se mantiene casi igual
# porque las perturbaciones son comunes). Es exacta como arco de ángulo en órbitas
# circulares; en elípticas es un desfase temporal constante, no una distancia.
#
# Funciones: along_track_layout
# =============================================================================

"""
    along_track_layout(chief, layout, ctx) -> Vector{Dict{String,Any}}

Genera `n` miembros en la órbita del jefe, desplazados en anomalía media
`spacing_km / SMA` cada uno respecto al anterior.

# Argumentos
- `chief::AbstractDict`: elementos keplerianos del jefe (de `build_orbit`).
- `layout::AbstractDict`: bloque `layout` con `n` (entero >= 1) y `spacing_km`.
- `ctx`: contexto del Builder.

# Devuelve
- `Vector` de `n` diccionarios de elementos; el primero es el jefe sin cambios.

# Excepciones
- `ErrorException`: si aparece una clave no admitida.
- `ErrorException`: si falta `n` o `spacing_km`, o `n` no es un entero >= 1.

# Ejemplo
```julia
chief = Dict("SMA" => 7000.0, "ECC" => 0.0, "INC" => 90.0, "RAAN" => 0.0, "AOP" => 0.0, "TA" => 0.0)
layout = Dict("type" => "along_track", "n" => 3, "spacing_km" => 10.0)
along_track_layout(chief, layout, ctx)      # TA = 0, 0.0819, 0.1637 grados (10 km de arco a 7000 km)
```
"""
function along_track_layout(chief::AbstractDict, layout::AbstractDict, ctx)
    check_keys(layout, ("type", "n", "spacing_km"); who = "along_track")

    nf = getreq(layout, "n")
    (isinteger(nf) && nf >= 1) || error("'along_track': 'n' must be an integer >= 1 (got $nf).")
    n = Int(nf)
    spacing = getreq(layout, "spacing_km")

    ecc  = chief["ECC"]
    dM   = rad2deg(spacing / chief["SMA"])                       # desfase en anomalía media entre vecinos (deg)
    M0   = true_to_mean_anomaly_deg(chief["TA"], ecc)

    members = Dict{String,Any}[]
    for k in 0:n-1
        el = Dict{String,Any}(chief)                             # copia: no se modifica el jefe
        # El jefe conserva su TA exacta; el resto se obtiene de M0 + k·ΔM
        k > 0 && (el["TA"] = mean_to_true_anomaly_deg(M0 + k * dM, ecc))
        push!(members, el)
    end
    return members
end

register_layout!("along_track", along_track_layout)
