# =============================================================================
# leo.jl
#
# Plugin de órbita: LEO genérica. Registra el tipo "LEO". Es keplerian con
# unidades más cómodas (alturas en vez de SMA) y una comprobación de apogeo.
#
# Dos formas excluyentes de dar el tamaño (alturas medidas sobre ctx.R):
#   a) circular:  {"type": "LEO", "altitude_km": 600, "INC": 51.6, "RAAN": 0, "TA": 0}
#        SMA = R + altitude_km,  ECC = 0
#   b) elíptica:  {"type": "LEO", "perigee_km": 500, "apogee_km": 900,
#                  "INC": 70, "RAAN": 0, "AOP": 90, "TA": 0}
#        SMA = R + (hp + ha)/2,  ECC = (ha - hp) / (ha + hp + 2R)
#
# Obligatorios: INC, RAAN y TA. AOP solo es obligatorio si la órbita es
# elíptica (en una circular no tiene significado físico y queda en 0).
#
# Funciones: leo_elements
# =============================================================================

"""
    leo_elements(p, ctx) -> Dict{String,Any}

Calcula los elementos keplerianos de una órbita LEO a partir de alturas,
inclinación, RAAN y posición en la órbita.

# Argumentos
- `p::AbstractDict`: bloque `orbit` con `altitude_km` o `perigee_km`+`apogee_km`,
  más `INC`, `RAAN`, `TA` y, si `ECC > 0`, `AOP`.
- `ctx`: contexto del Builder (usa `ctx.R`).

# Devuelve
- `Dict{String,Any}` con `SMA, ECC, INC, RAAN, AOP, TA` (km y grados).

# Excepciones
- `ErrorException`: si aparece una clave no admitida (por ejemplo `ECC` o `SMA`).
- `ErrorException`: si falta el tamaño de la órbita, `INC`, `RAAN`, `TA` o `AOP` con `ECC > 0`.
- `ErrorException`: si `apogee_km < perigee_km` o se mezclan las dos formas.

# Ejemplo
```julia
p = Dict("type" => "LEO", "altitude_km" => 600, "INC" => 51.6, "RAAN" => 0, "TA" => 0)
leo_elements(p, ctx)        # SMA = 6978.1363, ECC = 0, INC = 51.6, RAAN = 0, AOP = 0, TA = 0
```
"""
function leo_elements(p::AbstractDict, ctx)
    check_keys(p, ("type", "altitude_km", "perigee_km", "apogee_km", "INC", "RAAN", "AOP", "TA");
               who = "LEO")

    sma, ecc = sma_ecc_from_altitudes(p, ctx; who = "LEO")

    # Aviso si el apogeo sale de la región LEO (2000 km)
    apogee_alt = sma * (1 + ecc) - ctx.R
    apogee_alt > 2000 && @warn "LEO: apogee altitude = $(round(apogee_alt, digits=1)) km (> 2000 km)."

    return Dict{String,Any}(
        "SMA"  => sma,
        "ECC"  => ecc,
        "INC"  => getreq(p, "INC"),
        "RAAN" => getreq(p, "RAAN"),
        "AOP"  => get_aop(p, ecc),
        "TA"   => getreq(p, "TA")
    )
end

register_orbit!("LEO", leo_elements)