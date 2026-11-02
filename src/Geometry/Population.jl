# =============================================================================
# Population.jl
#
# Modelo estadístico de la población de objetos en órbita, usado para sortear la
# forma y las dimensiones de un objetivo cuando el usuario no las da.
# Se incluye desde Geometry.jl, por lo que comparte su espacio de nombres.
#
# Procedimiento (todo con el generador aleatorio del objetivo):
#   1. Clase de objeto (fragment, payload, rocket_body, mission_related), con
#      probabilidad proporcional al número de objetos de cada clase de tamaño ≥ L_min.
#   2. Forma dentro de la clase (SHAPE_WEIGHTS_BY_CLASS).
#   3. Tamaño de la clase: L_c (fragmentos y objetos de misión) o área media
#      (satélites y etapas), y proporciones de la forma; se escala la forma a ese tamaño.
#
# Unidades: metros y m². Parámetros: src/Parameters.jl (fuentes indicadas allí).
#
# Funciones:
#   fragment_cumulative   : N(>L_c) relativo de los fragmentos (ley de potencias por tramos)
#   sample_fragment_lc    : sorteo de L_c de un fragmento en [lo, hi]
#   sample_loguniform     : sorteo log-uniforme en [lo, hi]
#   sample_area_bins      : sorteo de un área a partir de un histograma por intervalos
#   class_weights         : número relativo de objetos de cada clase con L_c ≥ L_min
#   sample_class          : sorteo de la clase de objeto
#   sample_shape_type     : sorteo de la forma dentro de una clase
#   sample_size           : sorteo del tamaño (L_c o área) de una clase
#   sample_dims           : medidas (m) de una forma con el tamaño de una clase
# =============================================================================

# Clases de objeto que admite el modelo
const SHAPE_CLASSES = (:fragment, :payload, :rocket_body, :mission_related)

"""
    fragment_cumulative(L) -> Float64

Número acumulado relativo de fragmentos mayores que `L`, normalizado a 1 en
`L = 0.1 m`: `N(>L) ∝ L^-β` con el β de cada tramo de `FRAGMENT_LC_BREAKS_M`
(continuo en los límites de los tramos).

# Argumentos
- `L::Real`: longitud característica (m), dentro de los tramos.

# Devuelve
- `Float64`: N(>L) / N(>0.1 m).

# Excepciones
- `ErrorException`: si `L` está fuera de los tramos definidos.

# Ejemplo
```julia
fragment_cumulative(0.01)       # ≈ 22.9 (10^1.36): hay 23 veces más fragmentos de > 1 cm que de > 10 cm
```
"""
function fragment_cumulative(L::Real)
    b, β = FRAGMENT_LC_BREAKS_M, FRAGMENT_LC_BETAS
    (b[1] <= L <= b[end]) || error("fragment_cumulative: L = $L m outside [$(b[1]), $(b[end])] m.")

    # N(>b_i) en cada límite, empezando por N(>b_1) = 1 y encadenando los tramos
    Nb = ones(length(b))
    for i in eachindex(β)
        Nb[i + 1] = Nb[i] * (b[i + 1] / b[i])^(-β[i])
    end

    i = min(searchsortedlast(collect(b), L), length(β))   # tramo que contiene L
    N = Nb[i] * (L / b[i])^(-β[i])

    # Normalización en 0.1 m (límite del catálogo, donde están los recuentos por clase)
    i10 = min(searchsortedlast(collect(b), 0.1), length(β))
    return N / (Nb[i10] * (0.1 / b[i10])^(-β[i10]))
end

"""
    sample_fragment_lc(rng, lo, hi) -> Float64

Sortea la longitud característica de un fragmento con la ley de potencias por
tramos, truncada a `[lo, hi]`. Se elige el tramo con probabilidad proporcional a
su número de objetos y dentro de él se invierte la función de distribución:

`L = [a^-β - u (a^-β - b^-β)]^(-1/β)`, con `u ~ U(0, 1)` y `[a, b]` el tramo recortado.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `lo::Real`, `hi::Real`: límites del sorteo (m), dentro de `FRAGMENT_LC_BREAKS_M`.

# Devuelve
- `Float64`: L_c (m).

# Excepciones
- `ErrorException`: si `lo >= hi` o están fuera de los tramos.

# Ejemplo
```julia
L = sample_fragment_lc(MersenneTwister(1), 0.1, 1.0)    # casi siempre cerca de 0.1 m
```
"""
function sample_fragment_lc(rng::AbstractRNG, lo::Real, hi::Real)
    b, β = FRAGMENT_LC_BREAKS_M, FRAGMENT_LC_BETAS
    (b[1] <= lo < hi <= b[end]) || error("sample_fragment_lc: need $(b[1]) <= lo < hi <= $(b[end]) m (got $lo, $hi).")

    # Tramos recortados a [lo, hi] y su peso N(>a) - N(>b)
    segs = Tuple{Float64,Float64,Float64}[]
    w = Float64[]
    for i in eachindex(β)
        a, c = max(b[i], lo), min(b[i + 1], hi)
        a < c || continue
        push!(segs, (a, c, β[i]))
        push!(w, fragment_cumulative(a) - fragment_cumulative(c))
    end

    # Tramo elegido con probabilidad proporcional a su peso
    k = searchsortedfirst(cumsum(w) ./ sum(w), rand(rng))
    a, c, βk = segs[min(k, length(segs))]

    # Inversión de la función de distribución dentro del tramo
    u = rand(rng)
    return (a^(-βk) - u * (a^(-βk) - c^(-βk)))^(-1 / βk)
end

# Sorteo log-uniforme en [lo, hi]: tantos objetos por década en todo el rango
sample_loguniform(rng::AbstractRNG, lo::Real, hi::Real) = lo * (hi / lo)^rand(rng)

"""
    sample_area_bins(rng, bins) -> Float64

Sortea un área a partir de un histograma por intervalos: se elige el intervalo
con probabilidad proporcional a su número de objetos y dentro de él el área es
log-uniforme (los histogramas de ESA tienen el eje de áreas logarítmico).

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `bins`: colección de `(área mín, área máx, número de objetos)` en m².

# Devuelve
- `Float64`: área (m²).

# Excepciones
- Ninguna.
"""
function sample_area_bins(rng::AbstractRNG, bins)
    counts = [Float64(c) for (_, _, c) in bins]
    k = searchsortedfirst(cumsum(counts) ./ sum(counts), rand(rng))
    lo, hi, _ = bins[min(k, length(bins))]
    return sample_loguniform(rng, lo, hi)
end

"""
    class_weights(lc_min=SHAPE_LC_MIN_M) -> NamedTuple

Número relativo de objetos de cada clase con L_c ≥ `lc_min`. Los satélites, etapas
y objetos de misión son los catalogados de la tabla 3.1 de ESA (supuestos ≥ 10 cm);
los fragmentos se escalan con la ley de potencias desde los catalogados a 10 cm.

# Argumentos
- `lc_min::Real`: tamaño mínimo (m), dentro de `FRAGMENT_LC_BREAKS_M`.

# Devuelve
- `NamedTuple` con un peso por clase (`fragment`, `payload`, `rocket_body`, `mission_related`).

# Excepciones
- `ErrorException`: si `lc_min` está fuera de los tramos de los fragmentos.

# Ejemplo
```julia
class_weights(0.10)      # proporciones de la tabla 3.1 (≈ 34 % fragmentos)
class_weights(0.01)      # ≈ 92 % fragmentos
```
"""
function class_weights(lc_min::Real = SHAPE_LC_MIN_M)
    c = SHAPE_CLASS_COUNTS
    return (fragment        = c.fragment * fragment_cumulative(lc_min),
            payload         = Float64(c.payload),
            rocket_body     = Float64(c.rocket_body),
            mission_related = Float64(c.mission_related))
end

"""
    sample_class(rng; lc_min=SHAPE_LC_MIN_M) -> Symbol

Sortea la clase de objeto con probabilidad proporcional a `class_weights(lc_min)`.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `lc_min::Real`: tamaño mínimo (m).

# Devuelve
- `Symbol`: `:fragment`, `:payload`, `:rocket_body` o `:mission_related`.

# Excepciones
- Las de `class_weights`.
"""
function sample_class(rng::AbstractRNG; lc_min::Real = SHAPE_LC_MIN_M)
    w = class_weights(lc_min)
    p = cumsum(collect(Float64, values(w))) ./ sum(values(w))
    return keys(w)[min(searchsortedfirst(p, rand(rng)), length(p))]
end

"""
    sample_shape_type(rng, class) -> String

Sortea la forma de un objeto de la clase dada con los pesos de `SHAPE_WEIGHTS_BY_CLASS`.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `class::Symbol`: clase de objeto.

# Devuelve
- `String`: nombre de la forma (`"prism"`, `"cylinder"`...).

# Excepciones
- `ErrorException`: si la clase no existe.
"""
function sample_shape_type(rng::AbstractRNG, class::Symbol)
    class in SHAPE_CLASSES || error("Unknown object class '$class'. Available: $(join(SHAPE_CLASSES, ", ")).")
    w = getfield(SHAPE_WEIGHTS_BY_CLASS, class)
    p = cumsum(collect(Float64, values(w))) ./ sum(values(w))
    return String(keys(w)[min(searchsortedfirst(p, rand(rng)), length(p))])
end

"""
    sample_size(rng, class; lc_min=SHAPE_LC_MIN_M) -> NamedTuple

Sortea el tamaño de un objeto de la clase dada:
- `fragment`: L_c con la ley de potencias por tramos, en `[lc_min, L_max]`.
- `payload`, `rocket_body`: área media (m²) de los histogramas de ESA.
- `mission_related`: L_c log-uniforme en `MISSION_RELATED_LC_RANGE_M`.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `class::Symbol`: clase de objeto.
- `lc_min::Real`: tamaño mínimo (m).

# Devuelve
- `(kind, value)`: `kind = :lc` (valor en m) o `kind = :area` (valor en m², área media).

# Excepciones
- `ErrorException`: si la clase no existe o `lc_min` no deja rango para sortear.
"""
function sample_size(rng::AbstractRNG, class::Symbol; lc_min::Real = SHAPE_LC_MIN_M)
    if class === :fragment
        return (kind = :lc, value = sample_fragment_lc(rng, lc_min, FRAGMENT_LC_BREAKS_M[end]))
    elseif class === :payload
        return (kind = :area, value = sample_area_bins(rng, PAYLOAD_AREA_BINS_M2))
    elseif class === :rocket_body
        return (kind = :area, value = sample_area_bins(rng, ROCKET_BODY_AREA_BINS_M2))
    elseif class === :mission_related
        lo, hi = max(MISSION_RELATED_LC_RANGE_M[1], lc_min), MISSION_RELATED_LC_RANGE_M[2]
        lo < hi || error("SHAPE_LC_MIN_M = $lc_min m leaves no range for mission-related objects.")
        return (kind = :lc, value = sample_loguniform(rng, lo, hi))
    end
    error("Unknown object class '$class'. Available: $(join(SHAPE_CLASSES, ", ")).")
end

"""
    sample_dims(rng, type, class; lc_min=SHAPE_LC_MIN_M) -> Dict{String,Float64}

Medidas (m) de una forma con proporciones y tamaño aleatorios: se sortean las
proporciones (medidas a escala 1, función `unit` de la forma) y el tamaño de la
clase, y se escala la forma entera.

- Tamaño en L_c: factor `s = L_c / L_c(escala 1)`.
- Tamaño en área media: por el teorema de Cauchy, el área proyectada media de un
  cuerpo convexo es S/4 (S = superficie), así que `s = √(4 A / S(escala 1))`.

Si la L_c resultante es menor que `lc_min` (p. ej. un CubeSat muy pequeño), se repite
el sorteo: así todos los objetos sorteados cumplen L_c ≥ `lc_min`.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `type::AbstractString`: forma registrada.
- `class::Symbol`: clase de la que se toma el tamaño.
- `lc_min::Real`: tamaño mínimo (m).

# Devuelve
- `Dict{String,Float64}`: medidas de la forma (claves `*_m`).

# Excepciones
- `ErrorException`: si la forma no existe o no se logra L_c ≥ `lc_min` en 1000 intentos.

# Ejemplo
```julia
d = sample_dims(MersenneTwister(1), "cylinder", :rocket_body)    # Dict("radius_m" => ..., "height_m" => ...)
```
"""
function sample_dims(rng::AbstractRNG, type::AbstractString, class::Symbol; lc_min::Real = SHAPE_LC_MIN_M)
    entry = shape_entry(type)
    for _ in 1:1000
        unit = entry.unit(rng)                             # proporciones (escala 1)
        size = sample_size(rng, class; lc_min = lc_min)
        s = size.kind === :lc ? size.value / entry.lc(unit) :
                                sqrt(4 * size.value / entry.surface(unit))
        dims = Dict{String,Float64}(k => v * s for (k, v) in unit)
        entry.lc(dims) >= lc_min && return dims
    end
    error("sample_dims: could not draw a '$type' with L_c >= $lc_min m.")
end
