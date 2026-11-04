# =============================================================================
# Frames.jl
#
# Sistemas y direcciones de referencia que dependen de la órbita, y método
# TRIAD para construir una actitud a partir de dos direcciones. Se incluye
# desde Attitude.jl.
#
# LVLH (el mismo que Scenario/Relative_utils.jl, sistema de Hill):
#   x̂ = r/|r|          radial hacia fuera (zenit)
#   ẑ = h/|h|          normal al plano orbital, h = r × v
#   ŷ = ẑ × x̂          a lo largo de la traza (≈ velocidad en órbita casi circular)
#
# Direcciones de apuntado ("direction"), todas unitarias y en ECI:
#   "nadir" = -r̂   "zenith" = r̂   "velocity" = v̂   "orbit_normal" = ĥ   "sun" = ŝ
#   Un "-" delante invierte la dirección ("-sun", "-orbit_normal"...).
#   Un vector [x, y, z] es una dirección inercial fija en ECI.
#
# Funciones:
#   lvlh_to_eci     : DCM C_io (LVLH -> ECI) a partir de r y v
#   parse_direction : lee una "direction" de la spec
#   direction_vector: evalúa una "direction" en un instante
#   triad           : DCM C_ib que alinea dos ejes cuerpo con dos direcciones
# =============================================================================

"""
    lvlh_to_eci(r, v) -> SMatrix{3,3,Float64,9}

DCM que pasa componentes de LVLH a ECI (C_io). Sus columnas son los ejes
x̂ (radial), ŷ (a lo largo de la traza) y ẑ (normal orbital) expresados en ECI.

# Argumentos
- `r::AbstractVector`: posición ECI (km).
- `v::AbstractVector`: velocidad ECI (km/s).

# Devuelve
- `SMatrix{3,3,Float64,9}`: C_io, con `v_eci = C_io · v_lvlh`.

# Excepciones
- `ErrorException`: si `r` y `v` son paralelos (órbita degenerada, sin plano orbital).

# Ejemplo
```julia
C_io = lvlh_to_eci([7000, 0, 0], [0, 7.5, 0])     # identidad: x radial, y traza, z normal
```
"""
function lvlh_to_eci(r::AbstractVector, v::AbstractVector)
    h = cross(SVector{3,Float64}(r), SVector{3,Float64}(v))
    nh = norm(h)
    nh > 0 || error("LVLH frame undefined: position and velocity are parallel.")

    x̂ = SVector{3,Float64}(r) / norm(r)      # radial
    ẑ = h / nh                               # normal al plano orbital
    ŷ = cross(ẑ, x̂)                          # a lo largo de la traza
    return hcat(x̂, ŷ, ẑ)
end

# Nombres de dirección admitidos en "direction" (sin el signo)
const DIRECTION_NAMES = ("nadir", "zenith", "velocity", "orbit_normal", "sun")

"""
    parse_direction(x, who) -> NamedTuple

Lee una `direction` de la spec: un nombre de `DIRECTION_NAMES` (con `-` delante para
invertirlo) o un vector `[x, y, z]` fijo en ECI.

# Argumentos
- `x`: texto (p. ej. `"nadir"`, `"-sun"`) o vector de 3 números.
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.

# Devuelve
- `(kind::Symbol, sign::Float64, vec::SVector{3,Float64})`. `kind` es el nombre como
  `Symbol` o `:inertial`; `vec` solo se usa con `:inertial` (ya normalizado).

# Excepciones
- `ErrorException`: si el nombre no existe o el vector no es válido.

# Ejemplo
```julia
parse_direction("-sun", "deb")       # (kind = :sun, sign = -1.0, vec = [0, 0, 0])
parse_direction([0, 0, 2], "deb")    # (kind = :inertial, sign = 1.0, vec = [0, 0, 1])
```
"""
function parse_direction(x, who::AbstractString)
    if x isa AbstractString
        sgn  = startswith(x, "-") ? -1.0 : 1.0
        name = sgn < 0 ? x[2:end] : x        # nombre sin el signo
        name in DIRECTION_NAMES || error(
            "$who: unknown direction '$x'. Use one of $(join(DIRECTION_NAMES, ", ")) (optionally with '-') or an ECI vector [x, y, z].")
        return (kind = Symbol(name), sign = sgn, vec = zero(SVector{3,Float64}))
    end
    v = getunit(Dict{String,Any}("direction" => x), "direction", who)
    return (kind = :inertial, sign = 1.0, vec = v)
end

"""
    direction_vector(dir, r, v, s) -> SVector{3,Float64}

Evalúa una dirección de apuntado en un instante, en ECI.

# Argumentos
- `dir::NamedTuple`: dirección leída con `parse_direction`.
- `r::AbstractVector`: posición ECI del satélite (km).
- `v::AbstractVector`: velocidad ECI del satélite (km/s).
- `s::AbstractVector`: dirección unitaria del Sol en ECI.

# Devuelve
- `SVector{3,Float64}` unitario.

# Excepciones
- Ninguna.

# Ejemplo
```julia
direction_vector(parse_direction("nadir", "deb"), [7000, 0, 0], [0, 7.5, 0], [1, 0, 0])   # [-1, 0, 0]
```
"""
function direction_vector(dir, r::AbstractVector, v::AbstractVector, s::AbstractVector)
    r3 = SVector{3,Float64}(r)
    v3 = SVector{3,Float64}(v)

    u = dir.kind === :nadir        ? -r3 / norm(r3) :
        dir.kind === :zenith       ?  r3 / norm(r3) :
        dir.kind === :velocity     ?  v3 / norm(v3) :
        dir.kind === :orbit_normal ?  normalize(cross(r3, v3)) :
        dir.kind === :sun          ?  SVector{3,Float64}(s) / norm(s) :
                                      dir.vec                       # :inertial
    return dir.sign * u
end

"""
    triad(b1, b2, r1, r2) -> SMatrix{3,3,Float64,9}

Método TRIAD: DCM cuerpo -> ECI (C_ib) que lleva el eje cuerpo `b1` exactamente
a la dirección `r1` y deja `b2` lo más cerca posible de `r2`.

Con dos direcciones se forma una base ortonormal en cada sistema,
`t1 = b1`, `t2 = (b1 × b2)/|b1 × b2|`, `t3 = t1 × t2` (igual con `r1`, `r2`),
y entonces `C_ib = [t1r t2r t3r] · [t1b t2b t3b]ᵀ`. Una sola dirección no
basta: dejaría libre el giro alrededor de ella.

# Argumentos
- `b1, b2::AbstractVector`: ejes primario y secundario en ejes cuerpo (unitarios, no paralelos).
- `r1, r2::AbstractVector`: direcciones primaria y secundaria en ECI (unitarias).

# Devuelve
- `SMatrix{3,3,Float64,9}`: C_ib, con `C_ib · b1 = r1`.

# Excepciones
- `ErrorException`: si `r1` y `r2` son paralelas (la actitud queda indefinida).

# Ejemplo
```julia
# eje z del cuerpo al nadir y eje x hacia la velocidad
C_ib = triad([0, 0, 1], [1, 0, 0], -r/norm(r), v/norm(v))
```
"""
function triad(b1::AbstractVector, b2::AbstractVector, r1::AbstractVector, r2::AbstractVector)
    # Base del cuerpo
    t1b = SVector{3,Float64}(b1)
    t2b = normalize(cross(t1b, SVector{3,Float64}(b2)))
    t3b = cross(t1b, t2b)

    # Base en ECI; si las dos direcciones son paralelas no hay plano que fije el giro
    t1r = SVector{3,Float64}(r1)
    c   = cross(t1r, SVector{3,Float64}(r2))
    norm(c) > 1e-8 || error("Primary and secondary directions are parallel: the attitude is undefined.")
    t2r = c / norm(c)
    t3r = cross(t1r, t2r)

    return hcat(t1r, t2r, t3r) * hcat(t1b, t2b, t3b)'
end
