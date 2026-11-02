# =============================================================================
# fragment.jl
#
# Forma "fragment": fragmento irregular (envolvente convexa de 10–30 puntos aleatorios),
# el objeto más común de la población en órbita.
#
#   "shape": {"type": "fragment", "lc_m": 0.2}
#
# Ejes cuerpo: origen en el centroide; ejes principales, x la dimensión mayor (X del
# NASA SBM), y la siguiente y z la menor. La malla concreta depende de "seed".
# Sin medidas, L_c se sortea con la ley de potencias por tramos (clase fragment).
#
# Funciones: unit_fragment, lc_fragment, surface_fragment, mesh_fragment
# =============================================================================

# Medidas a escala 1: el fragmento solo tiene L_c
unit_fragment(rng::AbstractRNG) = Dict{String,Float64}("lc_m" => 1.0)

# L_c es la propia medida
lc_fragment(d) = d["lc_m"]

# Superficie aproximada de una esfera de diámetro L_c (solo para tamaños en área)
surface_fragment(d) = π * d["lc_m"]^2

# Envolvente convexa con la semilla del objetivo (medidas en km)
mesh_fragment(d, rng::AbstractRNG) = random_fragment(d["lc_m"]; rng = rng)

register_shape!("fragment"; dims = ("lc_m",), unit = unit_fragment, lc = lc_fragment,
                surface = surface_fragment, mesh = mesh_fragment, size_class = :fragment)
