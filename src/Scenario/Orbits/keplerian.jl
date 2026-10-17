# Órbita definida directamente por elementos keplerianos (km, grados).
# Reproduce exactamente lo que escribías a mano en los JSON.
#
#   {"type": "keplerian", "SMA": 8800, "ECC": 0.2, "INC": 90, "RAAN": 10, "AOP": 0, "TA": 357.5}

function keplerian_elements(p::AbstractDict, ctx)
    return Dict{String,Any}(
        "SMA"  => getreq(p, "SMA"),
        "ECC"  => getopt(p, "ECC", 0.0),
        "INC"  => getreq(p, "INC"),
        "RAAN" => getopt(p, "RAAN", 0.0),
        "AOP"  => getopt(p, "AOP", 0.0),
        "TA"   => getopt(p, "TA", 0.0)
    )
end

register_orbit!("keplerian", keplerian_elements)
