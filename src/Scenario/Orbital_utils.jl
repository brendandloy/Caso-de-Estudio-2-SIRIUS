# =============================================================================
# Orbital_utils.jl
#
# Funciones auxiliares compartidas por los plugins de órbita (Orbits/*.jl).
# Se incluye desde Scenario.jl, por lo que comparte su espacio de nombres
# (también ve getreq, getopt, ELEMENT_KEYS, etc.).
#
# Constantes:
#   J2_EARTH, TROPICAL_YEAR_DAYS
#
# Funciones:
#   check_keys              : rechaza claves no admitidas o calculadas.
#   sma_ecc_from_altitudes  : (altitud) o (perigeo, apogeo) -> (SMA, ECC).
#   get_aop                 : AOP obligatorio solo si la órbita no es circular.
#   parse_ltan              : "10:30" / 10.5 -> horas decimales.
#   ut_hours                : hora UT del día de una época, en horas.
#   era_deg                 : ángulo de rotación terrestre (ERA) del simulador.
#   sso_inclination_deg     : inclinación de una órbita heliosíncrona.
#   raan_from_ltan          : RAAN que corresponde a una LTAN dada.
#
# Importante: era_deg llama a Transformations.theta0_from_utc, la MISMA función
# que usa el simulador para rotar la Tierra. Así la longitud de un GEO y la
# LTAN de una SSO son coherentes con el simulador por construcción.
#
# Los ejemplos de los docstrings usan este contexto:
#   ctx = (R = 6378.1363, mu = 398600.4415, epoch = DateTime("2000-01-01T11:59:28"))
# =============================================================================

# J2 de la Tierra (EGM2008, sin normalizar); solo se usa para la inclinación SSO
const J2_EARTH = 1.08262668e-3

# Duración del año trópico en días; el nodo de una SSO debe seguir al Sol medio
const TROPICAL_YEAR_DAYS = 365.2421897

"""
    check_keys(p, allowed; who="orbit", computed=()) -> Nothing

Comprueba que el bloque `p` solo contiene claves admitidas por el tipo de órbita.

Por qué existe: sin ella, una errata como `"RAAm"` se ignoraría en silencio y el
plugin usaría otro valor sin avisar. También impide escribir a mano elementos que
el tipo calcula solo (p. ej. `INC` en una SSO).

# Argumentos
- `p::AbstractDict`: bloque `orbit` tal como viene de la spec.
- `allowed`: colección de claves admitidas (incluye `"type"`).
- `who::AbstractString`: nombre del tipo de órbita, solo para los mensajes de error.
- `computed`: claves que este tipo calcula por sí mismo; si aparecen, el error lo explica.

# Devuelve
- `nothing` si todas las claves son válidas.

# Excepciones
- `ErrorException`: si aparece una clave de `computed` (no se puede fijar a mano).
- `ErrorException`: si aparece cualquier otra clave que no esté en `allowed`.

# Ejemplo
```julia
p = Dict("type" => "GEO", "longitude_deg" => 0.0, "INC" => 0)
check_keys(p, ("type", "longitude_deg"); who = "GEO", computed = ("INC",))
# error: 'GEO': 'INC' is fixed/computed by this orbit type and cannot be set. ...
```
"""
function check_keys(p::AbstractDict, allowed; who::AbstractString = "orbit", computed = ())
    allowed_set  = Set(String.(collect(allowed)))
    computed_set = Set(String.(collect(computed)))

    for k in keys(p)
        ks = String(k)                  # las claves pueden ser String o Symbol según el parser
        ks in allowed_set && continue

        if ks in computed_set
            error("'$who': '$ks' is fixed/computed by this orbit type and cannot be set. " *
                  "Use type 'keplerian' if you need to define it manually.")
        end
        error("'$who': unknown parameter '$ks'. Allowed: " * join(sort!(collect(allowed_set)), ", "))
    end
    return nothing
end

"""
    sma_ecc_from_altitudes(p, ctx; who) -> Tuple{Float64,Float64}

Obtiene semieje mayor (km) y excentricidad a partir de alturas sobre `ctx.R`.
Admite dos formas excluyentes:
- `altitude_km`: órbita circular, `SMA = R + h`, `ECC = 0`.
- `perigee_km` + `apogee_km`: `SMA = R + (hp+ha)/2`, `ECC = (ha-hp)/(ha+hp+2R)`.

Por qué existe: `leo` y `sso` definen el tamaño de la órbita igual; así la
conversión está escrita una sola vez.

# Argumentos
- `p::AbstractDict`: bloque `orbit` de la spec.
- `ctx`: contexto del Builder (usa `ctx.R`, radio terrestre en km).
- `who::AbstractString`: nombre del tipo de órbita, para los mensajes de error.

# Devuelve
- `Tuple{Float64,Float64}`: `(sma_km, ecc)`.

# Excepciones
- `ErrorException`: si se mezclan ambas formas o falta alguna clave.
- `ErrorException`: si `apogee_km < perigee_km`.

# Ejemplo
```julia
sma_ecc_from_altitudes(Dict("altitude_km" => 600.0), ctx; who = "LEO")
# (6978.1363, 0.0)

sma_ecc_from_altitudes(Dict("perigee_km" => 614.8637, "apogee_km" => 628.8637), ctx; who = "LEO")
# (7000.0, 0.001)
```
"""
function sma_ecc_from_altitudes(p::AbstractDict, ctx; who::AbstractString)
    has_alt = haskey(p, "altitude_km")
    has_per = haskey(p, "perigee_km")
    has_apo = haskey(p, "apogee_km")

    if has_alt && (has_per || has_apo)
        error("'$who': use either 'altitude_km' (circular) or 'perigee_km' + 'apogee_km', not both.")
    elseif has_alt
        # Forma circular: la altura es constante
        return ctx.R + getreq(p, "altitude_km"), 0.0
    elseif has_per && has_apo
        hp = getreq(p, "perigee_km")
        ha = getreq(p, "apogee_km")
        ha >= hp || error("'$who': apogee_km must be >= perigee_km.")
        # Las alturas se miden sobre ctx.R, así que rp = R + hp y ra = R + ha
        return ctx.R + (hp + ha) / 2, (ha - hp) / (ha + hp + 2 * ctx.R)
    else
        error("'$who': provide 'altitude_km' or both 'perigee_km' and 'apogee_km'.")
    end
end

# AOP (deg): obligatorio si la órbita es elíptica; en una circular no tiene significado y vale 0
get_aop(p::AbstractDict, ecc::Real) = ecc > 0 ? getreq(p, "AOP") : getopt(p, "AOP", 0.0)

"""
    parse_ltan(x) -> Float64

Convierte una hora local del nodo ascendente (LTAN) a horas decimales.
Acepta un número (`10.5`) o un texto `"HH:MM"` / `"HH:MM:SS"`.

Por qué existe: en la spec es más cómodo escribir `"10:30"` que `10.5`, pero los
cálculos necesitan horas decimales.

# Argumentos
- `x`: LTAN como número de horas o como texto.

# Devuelve
- `Float64`: horas en el intervalo [0, 24).

# Excepciones
- `ErrorException`: si el texto no tiene formato válido o el valor queda fuera de [0, 24).

# Ejemplo
```julia
parse_ltan("10:30")       # 10.5
parse_ltan("10:30:36")    # 10.51
parse_ltan(14.25)         # 14.25
parse_ltan("25:00")       # error: LTAN must be in [0, 24) h (got 25.0).
```
"""
function parse_ltan(x)::Float64
    if x isa Real
        h = Float64(x)
    elseif x isa AbstractString
        parts = split(strip(x), ":")
        (1 <= length(parts) <= 3) || error("LTAN '$x': use 'HH:MM' or 'HH:MM:SS'.")

        nums = tryparse.(Float64, parts)
        any(isnothing, nums) && error("LTAN '$x': non-numeric field.")
        # Minutos y segundos deben estar en [0, 60)
        all(0 .<= nums[2:end] .< 60) || error("LTAN '$x': minutes and seconds must be in [0, 60).")

        h = nums[1] + (length(nums) >= 2 ? nums[2] / 60 : 0.0) + (length(nums) == 3 ? nums[3] / 3600 : 0.0)
    else
        error("LTAN must be a number of hours or a 'HH:MM' string (got $(typeof(x))).")
    end

    (0 <= h < 24) || error("LTAN must be in [0, 24) h (got $h).")
    return h
end

# Hora UT del día de la época, en horas decimales (p. ej. 11:59:28 -> 11.9911)
ut_hours(epoch::DateTime)::Float64 =
    hour(epoch) + minute(epoch) / 60 + (second(epoch) + millisecond(epoch) / 1000) / 3600

# ERA (deg) en la época dada; es el θ₀ del simulador (mismo datetime2julian que Propagator)
era_deg(epoch::DateTime)::Float64 = rad2deg(theta0_from_utc(datetime2julian(epoch)))

"""
    sso_inclination_deg(sma, ecc, ctx) -> Float64

Inclinación (deg) para que el plano orbital preceda a la velocidad media de
la Tierra alrededor del Sol. Usa la precesión nodal por J2:

    Ω̇ = -3/2 · J2 · (R/p)² · n · cos(i)   con   Ω̇ = 2π / año trópico

y despeja `cos(i) = -2Ω̇ / (3 J2 (R/p)² n)`, con `p = a(1-e²)` y `n = √(μ/a³)`.
Es una fórmula de primer orden en J2; la diferencia con elementos osculadores
es de décimas de grado y no afecta a la simulación.

# Argumentos
- `sma::Real`: semieje mayor en km.
- `ecc::Real`: excentricidad.
- `ctx`: contexto del Builder (usa `ctx.mu` y `ctx.R`).

# Devuelve
- `Float64`: inclinación en grados (> 90°, órbita retrógrada).

# Excepciones
- `ErrorException`: si `|cos(i)| > 1`, es decir, no existe SSO con esa altura y excentricidad
  (para órbita circular el límite está en torno a 5970 km de altitud).

# Ejemplo
```julia
sso_inclination_deg(ctx.R + 600, 0.0, ctx)      # ≈ 97.79
sso_inclination_deg(ctx.R + 700, 0.0, ctx)      # ≈ 98.19
sso_inclination_deg(ctx.R + 7000, 0.0, ctx)     # error: No sun-synchronous orbit exists ...
```
"""
function sso_inclination_deg(sma::Real, ecc::Real, ctx)::Float64
    n    = sqrt(ctx.mu / sma^3)                         # movimiento medio kepleriano (rad/s)
    p    = sma * (1 - ecc^2)                            # semilatus rectum (km)
    Ωdot = 2π / (TROPICAL_YEAR_DAYS * 86400.0)          # velocidad nodal requerida (rad/s)

    c = -2 * Ωdot / (3 * J2_EARTH * (ctx.R / p)^2 * n)  # cos(i)

    abs(c) <= 1 || error("No sun-synchronous orbit exists for SMA=$(round(sma, digits=1)) km, " *
                         "ECC=$ecc (cos(i) = $(round(c, digits=3)) is outside [-1, 1]).")
    return rad2deg(acos(c))
end

"""
    raan_from_ltan(ltan_h, epoch) -> Float64

RAAN (deg) cuyo nodo ascendente tiene hora solar media local `ltan_h` en `epoch`.
La longitud este del nodo es `λ = 15°·(LTAN − UT)`, y su ascensión recta es
`RAAN = ERA + λ`. La LTAN es hora solar MEDIA: el Sol verdadero (el que usa
SPICE para la iluminación) difiere de ella según la ecuación del tiempo, hasta
unos ±16 min en LTAN (≈ ±4° en RAAN) según la fecha.

Por qué existe: en una SSO el parámetro con sentido físico es la hora local a la
que cruza el ecuador (LTAN), no el RAAN; esta función hace la traducción.

# Argumentos
- `ltan_h::Real`: LTAN en horas decimales (ver `parse_ltan`).
- `epoch::DateTime`: época inicial del escenario (UTC).

# Devuelve
- `Float64`: RAAN en grados [0, 360).

# Excepciones
- Ninguna.

# Ejemplo
```julia
raan_from_ltan(10.5, DateTime("2000-01-01T11:59:28"))     # ≈ 257.96
```
"""
function raan_from_ltan(ltan_h::Real, epoch::DateTime)::Float64
    λ_node = 15.0 * (ltan_h - ut_hours(epoch))   # longitud este del nodo ascendente (deg)
    return mod(era_deg(epoch) + λ_node, 360.0)
end
