# =============================================================================
# random_leo.jl
#
# Plugin de órbita: LEO circular aleatoria. Registra el tipo "random_LEO".
#
#   {"type": "random_LEO"}                                       # todo aleatorio
#   {"type": "random_LEO", "altitude_km": [500, 800], "INC": 98}  # rango de altura, INC fija
#
# Cada parámetro admite tres formas:
#   - omitido          -> aleatorio en el rango por defecto
#   - [min, max]       -> aleatorio en ese rango
#   - número           -> valor fijo
#
# Rangos por defecto y distribución:
#   altitude_km  [400, 2000] km   uniforme (región LEO)
#   INC          [0, 180] deg     isótropa: cos(INC) uniforme, como si la normal al
#                                 plano orbital apuntara a una dirección al azar
#   RAAN         [0, 360] deg     uniforme
#   TA           [0, 360] deg     uniforme (en circular, argumento de latitud)
# ECC = 0 y AOP = 0 (circular); no se pueden escribir.
#
# Los números salen de ctx.rng, que el Builder crea con la "seed" de la spec: la
# misma semilla da siempre las mismas órbitas.
#
# Funciones: sample_param, random_leo_elements
# =============================================================================

"""
    sample_param(p, key, default, rng; isotropic=false) -> Float64

Devuelve el valor de un parámetro de una órbita aleatoria: fijo si es un número,
aleatorio en `[min, max]` si es un rango, o aleatorio en `default` si falta.

# Argumentos
- `p::AbstractDict`: bloque `orbit`.
- `key::AbstractString`: nombre del parámetro.
- `default::Tuple`: rango por defecto `(min, max)`.
- `rng`: generador aleatorio (`ctx.rng`).
- `isotropic::Bool`: si es `true`, el ángulo (deg) se muestrea con su coseno uniforme
  (se usa para la inclinación).

# Devuelve
- `Float64`: el valor elegido.

# Excepciones
- `ErrorException`: si el valor no es un número ni un rango `[min, max]` con `min <= max`.

# Ejemplo
```julia
rng = MersenneTwister(1)
sample_param(Dict("INC" => 98), "INC", (0, 180), rng)                    # 98.0
sample_param(Dict("altitude_km" => [500, 800]), "altitude_km", (400, 2000), rng)   # entre 500 y 800
```
"""
function sample_param(p::AbstractDict, key::AbstractString, default::Tuple, rng; isotropic::Bool = false)
    x = get(p, key, default)
    x isa Real && return Float64(x)                             # valor fijo

    (x isa Union{AbstractVector,Tuple} && length(x) == 2 && x[1] <= x[2]) ||
        error("'random_LEO': '$key' must be a number or a range [min, max] with min <= max.")
    lo, hi = Float64(x[1]), Float64(x[2])

    u = rand(rng)                                               # uniforme en [0, 1)
    if isotropic
        # cos(ángulo) uniforme entre cos(hi) y cos(lo): no se acumulan órbitas cerca de 0° y 180°
        return acosd(cosd(lo) + u * (cosd(hi) - cosd(lo)))
    end
    return lo + u * (hi - lo)
end

"""
    random_leo_elements(p, ctx) -> Dict{String,Any}

Genera los elementos de una órbita LEO circular aleatoria (ver cabecera del archivo
para rangos y distribuciones).

# Argumentos
- `p::AbstractDict`: bloque `orbit`; `altitude_km`, `INC`, `RAAN` y `TA` son opcionales.
- `ctx`: contexto del Builder (usa `ctx.R` y `ctx.rng`).

# Devuelve
- `Dict{String,Any}` con `SMA, ECC = 0, INC, RAAN, AOP = 0, TA` (km y grados).

# Excepciones
- `ErrorException`: si aparece `SMA`, `ECC`, `AOP` u otra clave no admitida.
- `ErrorException`: si un parámetro no es un número ni un rango válido.

# Ejemplo
```julia
random_leo_elements(Dict("type" => "random_LEO"), ctx)
# p. ej. SMA = 7512.3, ECC = 0, INC = 63.1, RAAN = 211.7, AOP = 0, TA = 40.2 (depende de la semilla)
```
"""
function random_leo_elements(p::AbstractDict, ctx)
    check_keys(p, ("type", "altitude_km", "INC", "RAAN", "TA");
               who = "random_LEO", computed = ("SMA", "ECC", "AOP"))

    # Orden fijo de muestreo: con la misma semilla salen siempre los mismos valores
    alt  = sample_param(p, "altitude_km", (400.0, 2000.0), ctx.rng)
    inc  = sample_param(p, "INC",  (0.0, 180.0), ctx.rng; isotropic = true)
    raan = sample_param(p, "RAAN", (0.0, 360.0), ctx.rng)
    ta   = sample_param(p, "TA",   (0.0, 360.0), ctx.rng)

    alt > 2000 && @warn "random_LEO: altitude = $(round(alt, digits=1)) km (> 2000 km, outside LEO)."

    return Dict{String,Any}(
        "SMA"  => ctx.R + alt,
        "ECC"  => 0.0,
        "INC"  => inc,
        "RAAN" => raan,
        "AOP"  => 0.0,
        "TA"   => ta
    )
end

register_orbit!("random_LEO", random_leo_elements)
