# =============================================================================
# sphere.jl
#
# Forma "sphere": esfera hueca de radio r (p. ej. esferas de calibración o tanques).
#
#   "shape": {"type": "sphere", "radius_m": 0.5}
#
# Ejes cuerpo: origen en el centro; z pasa por los polos de la malla (el brillo no
# depende de la orientación). L_c = 2r. Sin medidas, r se sortea con el tamaño de
# los objetos de misión (no hay clase propia de esferas en el informe de ESA).
#
# Funciones: unit_sphere, lc_sphere, surface_sphere, mesh_sphere
# =============================================================================

# Medidas a escala 1: la esfera no tiene proporciones que sortear
unit_sphere(rng::AbstractRNG) = Dict{String,Float64}("radius_m" => 1.0)

# L_c = (2r + 2r + 2r)/3 = 2r
lc_sphere(d) = 2 * d["radius_m"]

# Superficie 4πr²
surface_sphere(d) = 4π * d["radius_m"]^2

# Malla UV de 20 × 20 (medidas en km)
mesh_sphere(d, rng::AbstractRNG) = sphere_mesh(d["radius_m"], 20, 20)

register_shape!("sphere"; dims = ("radius_m",), unit = unit_sphere, lc = lc_sphere,
                surface = surface_sphere, mesh = mesh_sphere, size_class = :mission_related)
