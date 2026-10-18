# =============================================================================
# sso.jl
#
# Plugin de órbita: heliosíncrona (SSO). Registra el tipo "SSO".
#
# Entradas (todas en la spec):
#   - Tamaño: "altitude_km" (circular) o "perigee_km" + "apogee_km" (elíptica)
#   - "ltan" : hora local media del nodo ascendente, "10:30" o 10.5
#   - "TA"   : anomalía verdadera (argumento de latitud si la órbita es circular)
#   - "AOP"  : obligatorio solo si la órbita es elíptica
#
#   {"type": "SSO", "altitude_km": 600, "ltan": "10:30", "TA": 0}
#
# Se calculan (y NO se pueden escribir, da error):
#   INC  : de la precesión nodal por J2 (ver Orbital_utils.sso_inclination_deg)
#   RAAN : de la LTAN y la época (ver Orbital_utils.raan_from_ltan)
#
# Si quieres fijar INC y RAAN a mano, usa el tipo "keplerian".
#
# Funciones: sso_elements
# =============================================================================

"""
    sso_elements(p, ctx) -> Dict{String,Any}

Calcula los elementos keplerianos de una órbita heliosíncrona a partir de su
tamaño, la hora local del nodo ascendente y la posición en la órbita.

# Argumentos
- `p::AbstractDict`: bloque `orbit` con `altitude_km` o `perigee_km`+`apogee_km`,
  `ltan`, `TA` y, si `ECC > 0`, `AOP`.
- `ctx`: contexto del Builder (usa `ctx.mu`, `ctx.R` y `ctx.epoch`).

# Devuelve
- `Dict{String,Any}` con `SMA, ECC, INC, RAAN, AOP, TA` (km y grados).

# Excepciones
- `ErrorException`: si aparece `INC`, `RAAN` u otra clave no admitida.
- `ErrorException`: si falta `ltan`, `TA`, el tamaño de la órbita o `AOP` con `ECC > 0`.
- `ErrorException`: si no existe SSO con esa altura y excentricidad.
- `ErrorException`: si `ltan` no tiene un formato válido.

# Ejemplo
```julia
p = Dict("type" => "SSO", "altitude_km" => 600, "ltan" => "10:30", "TA" => 0)
sso_elements(p, ctx)
# SMA = 6978.1363, ECC = 0, INC ≈ 97.79, RAAN ≈ 257.96, AOP = 0, TA = 0 (época 2000-01-01T11:59:28)
```
"""
function sso_elements(p::AbstractDict, ctx)
    # INC y RAAN salen de la física del tipo: error si el usuario intenta fijarlos
    check_keys(p, ("type", "altitude_km", "perigee_km", "apogee_km", "ltan", "AOP", "TA");
               who = "SSO", computed = ("INC", "RAAN"))

    sma, ecc = sma_ecc_from_altitudes(p, ctx; who = "SSO")

    haskey(p, "ltan") || error("'SSO': missing required parameter 'ltan' (e.g. \"10:30\").")
    ltan_h = parse_ltan(p["ltan"])

    # Inclinación que hace que el nodo siga al Sol medio; falla si no existe SSO
    inc = sso_inclination_deg(sma, ecc, ctx)

    # RAAN que da la LTAN pedida en la época inicial del escenario
    raan = raan_from_ltan(ltan_h, ctx.epoch)

    return Dict{String,Any}(
        "SMA"  => sma,
        "ECC"  => ecc,
        "INC"  => inc,
        "RAAN" => raan,
        "AOP"  => get_aop(p, ecc),
        "TA"   => getreq(p, "TA")
    )
end

register_orbit!("SSO", sso_elements)
