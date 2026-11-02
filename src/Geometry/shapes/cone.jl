# =============================================================================
# cone.jl
#
# Forma "cone": superficie lateral de un cono, sin base (toberas, adaptadores).
#
#   "shape": {"type": "cone", "radius_m": 0.5, "height_m": 1.0}
#
# Ejes cuerpo: origen en el centro de la base; z es el eje, de la base (z = 0) al
# vértice (z = h). Sorteado: h/r uniforme en CONE_HR_RANGE. Sin medidas, el tamaño
# sale de los objetos de misión. L_c = (4r + h)/3.
#
# Funciones: unit_cone, lc_cone, surface_cone, mesh_cone
# =============================================================================

# Medidas a escala 1: r = 1 y h/r uniforme en CONE_HR_RANGE
unit_cone(rng::AbstractRNG) = Dict{String,Float64}("radius_m" => 1.0,
    "height_m" => CONE_HR_RANGE[1] + (CONE_HR_RANGE[2] - CONE_HR_RANGE[1]) * rand(rng))

# L_c = (2r + 2r + h)/3
lc_cone(d) = (4 * d["radius_m"] + d["height_m"]) / 3

# Superficie lateral πr√(r² + h²) (sin base)
surface_cone(d) = π * d["radius_m"] * hypot(d["radius_m"], d["height_m"])

# Superficie lateral con 50 divisiones (medidas en km)
mesh_cone(d, rng::AbstractRNG) = hollow_cone(d["radius_m"], d["height_m"])

register_shape!("cone"; dims = ("radius_m", "height_m"), unit = unit_cone, lc = lc_cone,
                surface = surface_cone, mesh = mesh_cone, size_class = :mission_related)
