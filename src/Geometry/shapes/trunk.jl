# =============================================================================
# trunk.jl
#
# Forma "trunk": superficie lateral de un tronco de cono, sin bases (adaptadores de
# carga, faldones entre etapas).
#
#   "shape": {"type": "trunk", "radius_m": 1.2, "top_radius_m": 0.6, "height_m": 1.0}
#
# Ejes cuerpo: origen en el centro de la base mayor; z es el eje, hacia la base menor
# (z = h). Sorteado: h/R en TRUNK_HR_RANGE y R/r en TRUNK_RADIUS_RATIO_RANGE. Sin
# medidas, el tamaño sale de los objetos de misión. L_c = (4R + h)/3.
#
# Funciones: check_trunk, unit_trunk, lc_trunk, surface_trunk, mesh_trunk
# =============================================================================

"""
    check_trunk(d, who)

Comprueba que el radio superior es menor que el de la base (si fueran iguales sería
un cilindro sin tapas).

# Argumentos
- `d`: medidas (m).
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si `top_radius_m >= radius_m`.
"""
function check_trunk(d, who::AbstractString)
    d["top_radius_m"] < d["radius_m"] ||
        error("$who: 'top_radius_m' must be smaller than 'radius_m' (base radius).")
    return nothing
end

"""
    unit_trunk(rng) -> Dict{String,Float64}

Medidas del tronco a escala 1: R = 1, r = 1/(R/r) con R/r uniforme en
`TRUNK_RADIUS_RATIO_RANGE` y h/R uniforme en `TRUNK_HR_RANGE`.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.

# Devuelve
- `Dict` con `radius_m`, `top_radius_m`, `height_m`.

# Excepciones
- Ninguna.
"""
function unit_trunk(rng::AbstractRNG)
    klo, khi = TRUNK_RADIUS_RATIO_RANGE
    hlo, hhi = TRUNK_HR_RANGE
    return Dict{String,Float64}("radius_m" => 1.0,
                                "top_radius_m" => 1 / (klo + (khi - klo) * rand(rng)),
                                "height_m" => hlo + (hhi - hlo) * rand(rng))
end

# L_c = (2R + 2R + h)/3
lc_trunk(d) = (4 * d["radius_m"] + d["height_m"]) / 3

# Superficie lateral π(R + r)√((R - r)² + h²) (sin bases)
surface_trunk(d) = π * (d["radius_m"] + d["top_radius_m"]) *
                   hypot(d["radius_m"] - d["top_radius_m"], d["height_m"])

# Superficie lateral con 50 divisiones (medidas en km)
mesh_trunk(d, rng::AbstractRNG) = hollow_cone(d["radius_m"], d["height_m"]; top_radius = d["top_radius_m"])

register_shape!("trunk"; dims = ("radius_m", "top_radius_m", "height_m"), unit = unit_trunk,
                lc = lc_trunk, surface = surface_trunk, mesh = mesh_trunk,
                size_class = :mission_related, check = check_trunk)
