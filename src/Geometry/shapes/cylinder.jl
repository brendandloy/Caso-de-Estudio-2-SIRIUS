# =============================================================================
# cylinder.jl
#
# Forma "cylinder": cilindro cerrado (etapas de cohete).
#
#   "shape": {"type": "cylinder", "radius_m": 1.5, "height_m": 8.0}
#
# Ejes cuerpo: origen en el centro; z es el eje de revolución.
# Sorteado: h/r uniforme en CYLINDER_HR_RANGE. Sin medidas, el tamaño sale del
# histograma de área de las etapas de ESA (clase rocket_body). L_c = (4r + h)/3.
#
# Funciones: unit_cylinder, lc_cylinder, surface_cylinder, mesh_cylinder
# =============================================================================

# Medidas a escala 1: r = 1 y h/r uniforme en CYLINDER_HR_RANGE
unit_cylinder(rng::AbstractRNG) = Dict{String,Float64}("radius_m" => 1.0,
    "height_m" => CYLINDER_HR_RANGE[1] + (CYLINDER_HR_RANGE[2] - CYLINDER_HR_RANGE[1]) * rand(rng))

# L_c = (2r + 2r + h)/3
lc_cylinder(d) = (4 * d["radius_m"] + d["height_m"]) / 3

# Superficie 2πr² + 2πrh (con las tapas)
surface_cylinder(d) = 2π * d["radius_m"] * (d["radius_m"] + d["height_m"])

# Malla con 20 divisiones en el contorno (medidas en km)
mesh_cylinder(d, rng::AbstractRNG) = cylinder_mesh(d["radius_m"], d["height_m"]; Nt = 20, Nz = 1)

register_shape!("cylinder"; dims = ("radius_m", "height_m"), unit = unit_cylinder, lc = lc_cylinder,
                surface = surface_cylinder, mesh = mesh_cylinder, size_class = :rocket_body)
