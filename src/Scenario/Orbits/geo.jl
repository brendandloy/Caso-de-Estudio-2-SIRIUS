# =============================================================================
# geo.jl
#
# Plugin de órbita: GEO (geoestacionaria). Registra el tipo "GEO".
#
# La única entrada es la longitud del punto subsatélite:
#
#   {"type": "GEO", "longitude_deg": -3.7}      # longitud este positiva
#
# Todo lo demás está fijado por la definición de GEO o se calcula:
#   SMA  : (μ / ω_earth²)^(1/3), periodo = día sidéreo del simulador
#   ECC  : 0          INC : 0         (definición de GEO)
#   RAAN : 0          AOP : 0         (indefinidos con e = 0, i = 0)
#   TA   : ERA(t0) + longitude_deg    (longitud verdadera inercial)
#
# Escribir cualquier elemento orbital en un bloque "GEO" da error. Para una
# órbita geosíncrona con e o i distintos de cero, usa "keplerian".
#
# Funciones: geo_elements
# =============================================================================

"""
    geo_elements(p, ctx) -> Dict{String,Any}

Calcula los elementos keplerianos de un satélite geoestacionario sobre la
longitud indicada. Con `e = 0` e `i = 0` solo importa `Ω+ω+ν`, que es la
longitud inercial del satélite: `ERA(t0) + longitud`. Se asigna entera a `TA`.

# Argumentos
- `p::AbstractDict`: bloque `orbit`; solo admite `type` y `longitude_deg` (grados este).
- `ctx`: contexto del Builder (usa `ctx.mu` y `ctx.epoch`).

# Devuelve
- `Dict{String,Any}` con `SMA, ECC, INC, RAAN, AOP, TA` (km y grados).

# Excepciones
- `ErrorException`: si aparece cualquier clave distinta de `type` y `longitude_deg`.
- `ErrorException`: si falta `longitude_deg` o está fuera de [-180, 360].

# Ejemplo
```julia
geo_elements(Dict("type" => "GEO", "longitude_deg" => -3.7), ctx)
# SMA ≈ 42164.17, ECC = 0, INC = 0, RAAN = 0, AOP = 0, TA ≈ 276.63 (época 2000-01-01T11:59:28)
```
"""
function geo_elements(p::AbstractDict, ctx)
    # Cualquier elemento orbital está fijado por el tipo: se rechaza aunque valga lo mismo
    check_keys(p, ("type", "longitude_deg"); who = "GEO", computed = ELEMENT_KEYS)

    lon = getreq(p, "longitude_deg")
    (-180 <= lon <= 360) || error("'GEO': longitude_deg must be in [-180, 360] (got $lon).")

    # Radio de la órbita con periodo igual a la rotación terrestre del simulador (ω_earth)
    sma = cbrt(ctx.mu / ω_earth^2)

    # Longitud inercial del satélite = ángulo de rotación terrestre + longitud geográfica
    ta = mod(era_deg(ctx.epoch) + lon, 360.0)

    return Dict{String,Any}(
        "SMA"  => sma,
        "ECC"  => 0.0,
        "INC"  => 0.0,
        "RAAN" => 0.0,
        "AOP"  => 0.0,
        "TA"   => ta
    )
end

register_orbit!("GEO", geo_elements)
