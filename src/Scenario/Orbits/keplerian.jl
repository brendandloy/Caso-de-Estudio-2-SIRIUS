# =============================================================================
# keplerian.jl
#
# Plugin de órbita: elementos keplerianos en bruto (km, grados). Registra el
# tipo "keplerian". Es la vía de escape: reproduce exactamente lo que se
# escribía a mano en los JSON y sirve para cualquier órbita sin tipo propio.
#
#   {"type": "keplerian", "SMA": 8800, "ECC": 0.2, "INC": 90, "RAAN": 10, "AOP": 0, "TA": 357.5}
#
# Los seis elementos son obligatorios: ningún valor por defecto es neutro,
# y un olvido en RAAN o TA colocaría el objeto en otro sitio sin avisar.
#
# Funciones: keplerian_elements
# =============================================================================

"""
    keplerian_elements(p, ctx) -> Dict{String,Any}

Devuelve los seis elementos keplerianos tal como se escriben en la spec.

# Argumentos
- `p::AbstractDict`: bloque `orbit` con `SMA, ECC, INC, RAAN, AOP, TA` (km y grados).
- `ctx`: contexto del Builder (no se usa; el contrato de plugin lo exige).

# Devuelve
- `Dict{String,Any}` con `SMA, ECC, INC, RAAN, AOP, TA` como `Float64`.

# Excepciones
- `ErrorException`: si falta alguno de los seis elementos.
- `ErrorException`: si aparece una clave que no sea `type` ni un elemento kepleriano.

# Ejemplo
```julia
p = Dict("type" => "keplerian", "SMA" => 8800, "ECC" => 0.2, "INC" => 90,
         "RAAN" => 10, "AOP" => 0, "TA" => 357.5)
keplerian_elements(p, ctx)        # los mismos valores, convertidos a Float64
```
"""
function keplerian_elements(p::AbstractDict, ctx)
    # Solo se admiten "type" y los seis elementos (detecta erratas como "RAAm")
    check_keys(p, ("type", ELEMENT_KEYS...); who = "keplerian")

    return Dict{String,Any}(
        "SMA"  => getreq(p, "SMA"),
        "ECC"  => getreq(p, "ECC"),
        "INC"  => getreq(p, "INC"),
        "RAAN" => getreq(p, "RAAN"),
        "AOP"  => getreq(p, "AOP"),
        "TA"   => getreq(p, "TA")
    )
end

register_orbit!("keplerian", keplerian_elements)
