module Geometry

# =============================================================================
# Módulo Geometry
#
# Forma de los objetivos (RSOs): lee y valida el bloque "shape" de cada objetivo,
# sortea forma y medidas cuando el usuario no las da y genera la malla de facetas
# (normales y áreas en ejes cuerpo) y las propiedades ópticas que usa la fotometría.
#
#   "shape": {"type": "cylinder", "radius_m": 1.5, "height_m": 8.0}   medidas concretas
#   "shape": {"type": "cylinder"}                                       medidas aleatorias
#   "shape": {"type": "random"}                                         forma y medidas aleatorias
#
# Scripts que contiene:
#   Geometry.jl   : este archivo; tipo ShapeSpec, registro de formas, lectura del
#                   bloque, sorteo (resolve_shape), ópticas y generate_debris
#   Meshes.jl     : mallas (prisma, esfera, cilindro, cono, fragmento) y normales
#   Population.jl : modelo de población (clases de ESA, tamaños) para los sorteos
#   shapes/       : un archivo por forma: sphere, prism, plate, cylinder, cone,
#                   trunk, fragment
#
# Medidas (m) y longitud característica L_c = (X + Y + Z)/3 (NASA SBM; X, Y, Z son las
# dimensiones según los ejes cuerpo):
#   sphere   : radius_m                          L_c = 2r
#   prism    : x_m, y_m, z_m                     L_c = (x + y + z)/3
#   plate    : x_m, y_m, thickness_m             L_c = (x + y + t)/3
#   cylinder : radius_m, height_m                L_c = (4r + h)/3
#   cone     : radius_m, height_m                L_c = (4r + h)/3
#   trunk    : radius_m, top_radius_m, height_m  L_c = (4R + h)/3
#   fragment : lc_m                              L_c = lc_m
#
# Ejes cuerpo de cada forma (los usa la actitud, que gira estos ejes respecto a ECI):
#   sphere   : origen en el centro; z pasa por los polos de la malla (el brillo no depende de los ejes)
#   prism    : origen en el centro; aristas x_m, y_m, z_m según x, y, z; z longitudinal (sorteado: z ≥ x ≥ y)
#   plate    : origen en el centro; x, y en el plano de la placa y z normal a ella (espesor)
#   cylinder : origen en el centro; z es el eje de revolución
#   cone     : origen en el centro de la base; z es el eje, de la base (z = 0) al vértice (z = h)
#   trunk    : origen en el centro de la base mayor; z es el eje, hacia la base menor (z = h)
#   fragment : origen en el centroide; ejes principales, x la dimensión mayor y z la menor
#
# Reglas del bloque "shape":
#   - Obligatorio en cada objetivo; "type" es una forma registrada o "random".
#   - Medidas: todas las de la forma o ninguna (ninguna = se sortean con el modelo de
#     población). Con medidas concretas no hay mínimo: si L_c < SHAPE_LC_MIN_M se avisa
#     de que probablemente no se verá en la imagen.
#   - "seed": semilla de la malla del fragmento y de las propiedades ópticas. El
#     generador de JSON pone una por miembro; sin ella se sortea una y se imprime.
#   - "class": solo informativa (la escribe el generador cuando sortea la forma).
#
# Contrato de una forma (para añadir una nueva basta con crear un archivo en shapes/):
#   register_shape!("nueva"; dims, unit, lc, surface, mesh, size_class, [check], [optics])
#     dims       : claves de las medidas, en metros
#     unit(rng)  : medidas a escala 1 con proporciones aleatorias (Dict)
#     lc(d)      : L_c (m) de unas medidas d
#     surface(d) : superficie (m²), para convertir un área media en tamaño
#     mesh(d_km, rng) -> (V, F) con las medidas en km
#     size_class : clase de la que se toma el tamaño cuando no hay medidas
#     check(d, who)  : validaciones propias (opcional)
#     optics(rng, d) : (C_abs, C_spec, C_diff) (opcional; por defecto las del cuerpo)
#
# Unidades: medidas del bloque en m; mallas en km y áreas en km² (las de la fotometría).
# Dependencias: Parameters (parámetros maestros de la forma) y Attitude (lectores de
# bloques todict, check_keys, getnum), que deben estar cargados antes.
# =============================================================================

export ShapeSpec, Debris, register_shape!, shape_block, parse_shape, parse_member_shape
export resolve_shape, char_length, generate_debris

using ..Parameters
using ..Attitude: todict, check_keys, getnum      # mismos lectores de bloques que la actitud

using LinearAlgebra
using Random
using Statistics

# Diccionario que recibe la fotometría (normales, áreas, coeficientes ópticos...)
const Debris = Dict{Symbol,Any}

"""
    ShapeSpec(type, dims_m, seed)

Forma de un objetivo ya validada, con todas sus medidas.

# Campos
- `type::String`: forma registrada (`"prism"`, `"cylinder"`...).
- `dims_m::Dict{String,Float64}`: medidas en metros (claves de la forma, p. ej. `"radius_m"`).
- `seed::Int`: semilla de la malla (fragmento) y de las propiedades ópticas.
"""
struct ShapeSpec
    type::String
    dims_m::Dict{String,Float64}
    seed::Int
end

# -----------------------------------------------------------------------------
# Propiedades ópticas
# -----------------------------------------------------------------------------

# (μ_abs, σ_abs, μ_spec_aux, σ_spec_aux) de cada material (valores heredados del CE2 de Ortiz)
const OPTICS_BODY  = (0.62, 0.12, 0.7, 0.1)        # cuerpo (MLI, aluminio...)
const OPTICS_TOP   = (0.55, 0.2, 0.7, 0.1)         # tapas ±z del prisma
const OPTICS_PANEL = (0.75, 0.1, 0.8, 0.05)        # paneles y placas

"""
    rand_trunc_normal(rng, μ, σ; lower=0.0, upper=1.0) -> Float64

Normal N(μ, σ) truncada a `[lower, upper)` por rechazo.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `μ`, `σ`: media y desviación típica.
- `lower`, `upper`: límites del intervalo.

# Devuelve
- `Float64` dentro de `[lower, upper)`.

# Excepciones
- Ninguna (si el intervalo es muy improbable, el bucle tarda más).
"""
function rand_trunc_normal(rng::AbstractRNG, μ, σ; lower = 0.0, upper = 1.0)
    x = μ + σ * randn(rng)
    while (x < lower) || (x >= upper)
        x = μ + σ * randn(rng)
    end
    return x
end

"""
    optical_coeffs(rng, params) -> (C_abs, C_spec, C_diff)

Coeficientes de absorción, reflexión especular y difusa de un material (CE2 de Ortiz,
§4.4.3): `C_abs ~ N(μ_abs, σ_abs)` y `C_spec,aux ~ N(μ_spec, σ_spec)` truncadas a [0, 1),
`C_spec = C_spec,aux (1 - C_abs)` y `C_diff = 1 - C_abs - C_spec`.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `params`: `(μ_abs, σ_abs, μ_spec_aux, σ_spec_aux)`.

# Devuelve
- `(C_abs, C_spec, C_diff)`, que suman 1.

# Excepciones
- Ninguna.
"""
function optical_coeffs(rng::AbstractRNG, params)
    μa, σa, μs, σs = params
    Cabs = rand_trunc_normal(rng, μa, σa)
    Cspec = rand_trunc_normal(rng, μs, σs) * (1 - Cabs)
    Cdiff = max(1 - Cabs - Cspec, 0.0)
    return Cabs, Cspec, Cdiff
end

# Ópticas por defecto de una forma: material del cuerpo en todas las caras
body_optics(rng::AbstractRNG, d) = optical_coeffs(rng, OPTICS_BODY)

# -----------------------------------------------------------------------------
# Mallas y modelo de población (los usa el registro y el sorteo)
# -----------------------------------------------------------------------------

include(joinpath(@__DIR__, "Meshes.jl"))
include(joinpath(@__DIR__, "Population.jl"))

# -----------------------------------------------------------------------------
# Registro de formas
# -----------------------------------------------------------------------------

# Nombre de la forma -> funciones y datos que la definen (ver el contrato arriba)
const SHAPES = Dict{String,NamedTuple}()

# Comprobación propia por defecto (ninguna)
_no_check(d, who) = nothing

"""
    register_shape!(name; dims, unit, lc, surface, mesh, size_class, check, optics)

Registra una forma. Lo llama cada archivo de `shapes/` al final.

# Argumentos
- `name::AbstractString`: nombre de la forma en el bloque `shape` (`"type"`).
- Palabras clave: las del contrato de una forma (cabecera del módulo).

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si `size_class` no es una clase del modelo de población.

# Ejemplo
```julia
register_shape!("sphere"; dims = ("radius_m",), unit = rng -> Dict("radius_m" => 1.0),
                lc = d -> 2d["radius_m"], surface = d -> 4π * d["radius_m"]^2,
                mesh = (d, rng) -> sphere_mesh(d["radius_m"], 20, 20), size_class = :mission_related)
```
"""
function register_shape!(name::AbstractString; dims, unit, lc, surface, mesh, size_class::Symbol,
                         check = _no_check, optics = body_optics)
    size_class in SHAPE_CLASSES || error("register_shape!: unknown size_class '$size_class'.")
    SHAPES[String(name)] = (dims = Tuple(String.(dims)), unit = unit, lc = lc, surface = surface,
                            mesh = mesh, size_class = size_class, check = check, optics = optics)
    return nothing
end

# Datos registrados de una forma; error con la lista de formas si no existe
shape_entry(type::AbstractString) = haskey(SHAPES, type) ? SHAPES[type] :
    error("Unknown shape type '$type'. Available: $(join(sort!(collect(keys(SHAPES))), ", ")), random")

# -----------------------------------------------------------------------------
# Lectura del bloque "shape"
# -----------------------------------------------------------------------------

"""
    shape_block(entry, name) -> Dict{String,Any}

Devuelve el bloque `shape` de un objetivo (del JSON del simulador o de un grupo de
la spec) tras comprobar que existe, es un bloque `{...}` y tiene `type`.

# Argumentos
- `entry::AbstractDict`: entrada del objetivo (o grupo de objetivos de la spec).
- `name::AbstractString`: nombre del objetivo, para los mensajes de error.

# Devuelve
- `Dict{String,Any}`: copia del bloque `shape` (sin validar sus parámetros).

# Excepciones
- `ErrorException`: si falta el bloque `shape`, no es un bloque `{...}` o no tiene `type`.

# Ejemplo
```julia
shape_block(Dict("shape" => Dict("type" => "random")), "deb")     # Dict("type" => "random")
```
"""
function shape_block(entry, name::AbstractString)
    e = todict(entry)
    haskey(e, "shape") || error(
        "Objective '$name': missing 'shape' block (e.g. \"shape\": {\"type\": \"random\"} or " *
        "\"shape\": {\"type\": \"cylinder\", \"radius_m\": 1.5, \"height_m\": 8.0}).")
    e["shape"] isa Dict{String,Any} || error("Objective '$name': 'shape' must be a block {...}.")
    haskey(e["shape"], "type") || error("Objective '$name': 'shape' block without 'type'.")
    return e["shape"]
end

"""
    getseed(b, who) -> Int

Lee la `seed` de un bloque.

# Argumentos
- `b::AbstractDict`: bloque con la clave `seed`.
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.

# Devuelve
- `Int` >= 0.

# Excepciones
- `ErrorException`: si `seed` no es un entero >= 0.
"""
function getseed(b::AbstractDict, who::AbstractString)
    s = b["seed"]
    (s isa Real && !(s isa Bool) && isinteger(s) && s >= 0) ||
        error("$who: 'seed' must be an integer >= 0 (got $(repr(s))).")
    return Int(s)
end

# Medidas que trae el bloque, de entre las de la forma
given_dims(b::AbstractDict, entry) = [k for k in entry.dims if haskey(b, k)]

"""
    parse_shape(block, name) -> ShapeSpec

Lee y valida un bloque `shape` con todas sus medidas (el que escribe el generador o
uno escrito a mano con medidas concretas). Si L_c es menor que `SHAPE_LC_MIN_M`, avisa
de que el objetivo probablemente no se verá en la imagen, pero lo acepta.

# Argumentos
- `block::AbstractDict`: bloque `shape`.
- `name::AbstractString`: nombre del objetivo, para los mensajes.

# Devuelve
- `ShapeSpec`.

# Excepciones
- `ErrorException`: si `type` no es una forma registrada (o es `random`, que hay que
  resolver antes con `resolve_shape`), hay claves desconocidas, falta alguna medida,
  alguna no es positiva, falla una comprobación propia de la forma o `seed` no es válida.

# Ejemplo
```julia
parse_shape(Dict("type" => "cylinder", "radius_m" => 1.5, "height_m" => 8.0, "seed" => 1), "deb")
```
"""
function parse_shape(block, name::AbstractString)
    b = todict(block)
    type = String(b["type"])
    who = "Objective '$name' (shape '$type')"
    type == "random" && error("$who: a 'random' shape must be resolved first (resolve_shape).")
    entry = shape_entry(type)

    # Claves admitidas y medidas completas
    check_keys(b, ("type", "class", "seed", entry.dims...), who)
    missing_dims = [k for k in entry.dims if !haskey(b, k)]
    isempty(missing_dims) || error("$who: missing $(join(missing_dims, ", ")) " *
        "(give all of $(join(entry.dims, ", ")) or none to draw them at random).")

    # Medidas positivas y comprobaciones propias de la forma
    dims = Dict{String,Float64}()
    for k in entry.dims
        v = getnum(b, k, who)
        v > 0 || error("$who: '$k' must be > 0 (got $v).")
        dims[k] = v
    end
    entry.check(dims, who)

    # Clase informativa: tiene que ser una de las del modelo de población
    haskey(b, "class") && !(Symbol(b["class"]) in SHAPE_CLASSES) &&
        error("$who: unknown 'class' '$(b["class"])'. Available: $(join(SHAPE_CLASSES, ", ")).")

    # Semilla de la malla y de las ópticas: la del bloque o una nueva que se imprime
    if haskey(b, "seed")
        seed = getseed(b, who)
    else
        seed = Int(rand(RandomDevice(), UInt32))
        println("Objective '$name': shape without 'seed'; using seed = $seed ",
                "(add \"seed\": $seed to its 'shape' block to repeat it).")
    end

    # Aviso (no error) si es más pequeño que el mínimo de los sorteos
    lc = entry.lc(dims)
    lc < SHAPE_LC_MIN_M * (1 - 1e-6) && @warn "$who: L_c = $(round(lc; sigdigits = 3)) m is below " *
        "SHAPE_LC_MIN_M = $SHAPE_LC_MIN_M m; such a small object will probably be too faint " *
        "to appear in the images unless it is very close to the camera."

    return ShapeSpec(type, dims, seed)
end

"""
    resolve_shape(block, name; seed) -> Dict{String,Any}

Convierte un bloque `shape` en otro con todas las medidas, sorteando lo que falte con
el modelo de población y la semilla dada:
- `type = "random"`: se sortean la clase, la forma y las medidas (se añade `class`).
- forma sin medidas: se sortean las medidas con el tamaño de la clase de esa forma.
- forma con todas las medidas: se copian tal cual.
En todos los casos se añade `seed`. El resultado se valida con `parse_shape`.

# Argumentos
- `block::AbstractDict`: bloque `shape` (de la spec o del JSON).
- `name::AbstractString`: nombre del objetivo, para los mensajes.
- `seed::Integer`: semilla del sorteo y de la malla/ópticas del objetivo.

# Devuelve
- `Dict{String,Any}`: bloque listo para el JSON del simulador.

# Excepciones
- `ErrorException`: forma desconocida, medidas incompletas, claves desconocidas o
  cualquier error de `parse_shape`.

# Ejemplo
```julia
resolve_shape(Dict("type" => "random"), "deb1"; seed = 42)
# Dict("type" => "cylinder", "class" => "rocket_body", "radius_m" => ..., "height_m" => ..., "seed" => 42)
```
"""
function resolve_shape(block, name::AbstractString; seed::Integer)
    b = todict(block)
    type = String(b["type"])
    who = "Objective '$name' (shape '$type')"
    rng = MersenneTwister(seed)
    out = Dict{String,Any}()

    if type == "random"
        # Forma y medidas aleatorias según la población real
        check_keys(b, ("type", "seed"), who)
        class = sample_class(rng)
        type = sample_shape_type(rng, class)
        merge!(out, round_dims(sample_dims(rng, type, class)))
        out["class"] = String(class)
    else
        entry = shape_entry(type)
        if isempty(given_dims(b, entry))
            # Forma elegida, medidas aleatorias con el tamaño de su clase
            check_keys(b, ("type", "seed"), who)
            merge!(out, round_dims(sample_dims(rng, type, entry.size_class)))
        else
            # Medidas concretas: se copian (parse_shape avisa si falta alguna)
            for (k, v) in b
                k in ("type", "seed") || (out[k] = v)
            end
        end
    end

    out["type"] = type
    out["seed"] = Int(seed)
    parse_shape(out, name)              # valida (y avisa si L_c < SHAPE_LC_MIN_M)
    return out
end

# Medidas sorteadas redondeadas a 6 cifras significativas, para que el JSON sea legible
round_dims(d) = Dict{String,Any}(k => round(v; sigdigits = 6) for (k, v) in d)

"""
    parse_member_shape(entry, name) -> ShapeSpec

Lee la forma de un objetivo del JSON del simulador. Si el bloque aún tiene algo por
sortear (`random` o sin medidas, p. ej. un JSON escrito a mano), lo resuelve con su
`seed` o con una nueva que se imprime.

# Argumentos
- `entry::AbstractDict`: entrada del objetivo en la sección `objectives`.
- `name::AbstractString`: nombre del objetivo, para los mensajes.

# Devuelve
- `ShapeSpec`.

# Excepciones
- Las de `shape_block`, `resolve_shape` y `parse_shape`.

# Ejemplo
```julia
entry = Dict("SMA" => 7000, "shape" => Dict("type" => "random", "seed" => 3))
parse_member_shape(entry, "deb")       # ShapeSpec con forma y medidas sorteadas
```
"""
function parse_member_shape(entry, name::AbstractString)
    b = shape_block(entry, name)
    type = String(b["type"])
    needs_draw = type == "random" || isempty(given_dims(b, shape_entry(type)))

    if needs_draw
        who = "Objective '$name' (shape '$type')"
        seed = haskey(b, "seed") ? getseed(b, who) : Int(rand(RandomDevice(), UInt32))
        haskey(b, "seed") || println("Objective '$name': shape to draw without 'seed'; using seed = $seed ",
                                     "(add \"seed\": $seed to its 'shape' block to repeat it).")
        b = resolve_shape(b, name; seed = seed)
    end
    return parse_shape(b, name)
end

# L_c (m) de una forma ya validada
char_length(s::ShapeSpec) = shape_entry(s.type).lc(s.dims_m)

# -----------------------------------------------------------------------------
# Malla y modelo fotométrico de un objetivo
# -----------------------------------------------------------------------------

"""
    generate_debris(spec::ShapeSpec) -> Debris

Genera el modelo de un objetivo para la fotometría: malla en ejes cuerpo, normales
y áreas de las facetas, coeficientes ópticos y exponente especular `alpha`. Todo lo
aleatorio (malla del fragmento, ópticas, `alpha`) usa `spec.seed`, así que la misma
spec da siempre el mismo objeto.

# Argumentos
- `spec::ShapeSpec`: forma validada (`parse_shape` o `parse_member_shape`).

# Devuelve
- `Debris` (`Dict{Symbol,Any}`) con:
  `:shape_name` (Symbol), `:dims_m`, `:char_length_m`, `:seed`,
  `:normals` (3 × N, ejes cuerpo), `:areas` (km²), `:total_area` (km²),
  `:C_abs`, `:C_spec`, `:C_diff` y `:alpha`.

# Excepciones
- Las de la malla de la forma (p. ej. `calculateNormals` si fuera degenerada).

# Ejemplo
```julia
spec = parse_shape(Dict("type" => "sphere", "radius_m" => 0.5, "seed" => 1), "deb")
d = generate_debris(spec)
d[:total_area]          # ≈ 4π (0.5e-3)² km²
```
"""
function generate_debris(spec::ShapeSpec)
    entry = shape_entry(spec.type)
    rng = MersenneTwister(spec.seed)

    # Malla en km (unidades de la fotometría) y facetas en ejes cuerpo
    dims_km = Dict{String,Float64}(k => 1e-3 * v for (k, v) in spec.dims_m)
    V, F = entry.mesh(dims_km, rng)
    normals, areas, total_area = calculateNormals(V, F)

    # Propiedades ópticas y exponente especular (Phong), con la misma semilla
    Cabs, Cspec, Cdiff = entry.optics(rng, spec.dims_m)
    alpha = max(0.1, 1.2 + 0.1 * randn(rng))

    debris = Debris()
    debris[:shape_name]    = Symbol(spec.type)
    debris[:dims_m]        = copy(spec.dims_m)
    debris[:char_length_m] = entry.lc(spec.dims_m)
    debris[:seed]          = spec.seed
    debris[:normals]       = normals
    debris[:areas]         = areas
    debris[:total_area]    = total_area
    debris[:C_abs]         = Cabs
    debris[:C_spec]        = Cspec
    debris[:C_diff]        = Cdiff
    debris[:alpha]         = alpha
    return debris
end

# -----------------------------------------------------------------------------
# Autocarga de formas
# -----------------------------------------------------------------------------

# Se incluyen todos los .jl de shapes/, en orden alfabético
let dir = joinpath(@__DIR__, "shapes")
    # Si el módulo se precompila, esto permite detectar archivos nuevos; si falla, no importa
    try
        include_dependency(dir)
    catch
    end
    for f in sort!(filter(endswith(".jl"), readdir(dir)))
        include(joinpath(dir, f))
    end
end

end # module Geometry
