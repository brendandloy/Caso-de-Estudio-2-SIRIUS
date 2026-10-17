# Órbita LEO genérica. Dos formas de definir el tamaño:
#
#   a) altitud media:      {"type": "leo", "altitude_km": 600, "ECC": 0.001, "INC": 51.6}
#        SMA = R + altitude_km  (ECC opcional, por defecto 0)
#   b) perigeo y apogeo:   {"type": "leo", "perigee_km": 500, "apogee_km": 900, "INC": 70}
#        SMA = R + (hp + ha)/2 ,  ECC = (ha - hp) / (ha + hp + 2R)
#
# RAAN, AOP y TA son opcionales (por defecto 0). Alturas medidas sobre R = ctx.R.

function leo_elements(p::AbstractDict, ctx)
    Re = ctx.R

    if haskey(p, "altitude_km")
        sma = Re + getreq(p, "altitude_km")
        ecc = getopt(p, "ECC", 0.0)
    elseif haskey(p, "perigee_km") && haskey(p, "apogee_km")
        hp = getreq(p, "perigee_km")
        ha = getreq(p, "apogee_km")
        ha >= hp || error("LEO: apogee_km must be >= perigee_km.")
        sma = Re + (hp + ha) / 2
        ecc = (ha - hp) / (ha + hp + 2Re)
    else
        error("LEO: provide 'altitude_km' or both 'perigee_km' and 'apogee_km'.")
    end

    apogee_alt = sma * (1 + ecc) - Re
    apogee_alt > 2000 && @warn "LEO: apogee altitude = $(round(apogee_alt, digits=1)) km (> 2000 km)."

    return Dict{String,Any}(
        "SMA"  => sma,
        "ECC"  => ecc,
        "INC"  => getreq(p, "INC"),
        "RAAN" => getopt(p, "RAAN", 0.0),
        "AOP"  => getopt(p, "AOP", 0.0),
        "TA"   => getopt(p, "TA", 0.0)
    )
end

register_orbit!("leo", leo_elements)
