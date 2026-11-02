# =============================================================================
# prism.jl
#
# Forma "prism": prisma rectangular (cuerpo de un satélite, sin paneles).
#
#   "shape": {"type": "prism", "x_m": 0.1, "y_m": 0.1, "z_m": 0.34}
#
# Ejes cuerpo: origen en el centro; aristas x_m, y_m, z_m según x, y, z. Convención:
# z es la dirección longitudinal (arista mayor) y las caras ±z son las tapas, con su
# propio material. Sorteado: z ≥ x ≥ y, con las aristas entre 1 y 3 veces la menor
# (PRISM_RATIO_RANGE). Sin medidas, el tamaño
# sale del histograma de área de los satélites de ESA (clase payload).
# L_c = (x + y + z)/3.
#
# Funciones: unit_prism, lc_prism, surface_prism, mesh_prism, optics_prism
# =============================================================================

"""
    unit_prism(rng) -> Dict{String,Float64}

Medidas del prisma a escala 1: arista menor 1 y las otras dos uniformes en
`PRISM_RATIO_RANGE`, ordenadas para que z ≥ x ≥ y (z es la dirección longitudinal).

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.

# Devuelve
- `Dict` con `x_m`, `y_m`, `z_m`.

# Excepciones
- Ninguna.
"""
function unit_prism(rng::AbstractRNG)
    lo, hi = PRISM_RATIO_RANGE
    a = sort([1.0, lo + (hi - lo) * rand(rng), lo + (hi - lo) * rand(rng)]; rev = true)
    # Arista mayor según z (longitudinal), la siguiente según x y la menor según y
    return Dict{String,Float64}("z_m" => a[1], "x_m" => a[2], "y_m" => a[3])
end

# L_c = (x + y + z)/3
lc_prism(d) = (d["x_m"] + d["y_m"] + d["z_m"]) / 3

# Superficie 2(xy + yz + xz)
surface_prism(d) = 2 * (d["x_m"] * d["y_m"] + d["y_m"] * d["z_m"] + d["x_m"] * d["z_m"])

# Malla de 12 triángulos (medidas en km)
mesh_prism(d, rng::AbstractRNG) = cube_vertices_faces(d["x_m"], d["y_m"], d["z_m"])

"""
    optics_prism(rng, d) -> (C_abs, C_spec, C_diff)

Coeficientes ópticos del prisma: media ponderada por área de los de las tapas ±z
(A_top = 2xy) y los del cuerpo (A_body = 2(xz + yz)), como en el CE2 de Ortiz.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.
- `d`: medidas (m).

# Devuelve
- `(C_abs, C_spec, C_diff)`.

# Excepciones
- Ninguna.
"""
function optics_prism(rng::AbstractRNG, d)
    _, Cs_b, Cd_b = optical_coeffs(rng, OPTICS_BODY)        # cuerpo
    _, Cs_t, Cd_t = optical_coeffs(rng, OPTICS_TOP)         # tapas ±z
    A_top  = 2 * d["x_m"] * d["y_m"]
    A_body = 2 * (d["x_m"] * d["z_m"] + d["y_m"] * d["z_m"])
    Cspec = (Cs_t * A_top + Cs_b * A_body) / (A_top + A_body)
    Cdiff = (Cd_t * A_top + Cd_b * A_body) / (A_top + A_body)
    return 1 - Cspec - Cdiff, Cspec, Cdiff
end

register_shape!("prism"; dims = ("x_m", "y_m", "z_m"), unit = unit_prism, lc = lc_prism,
                surface = surface_prism, mesh = mesh_prism, size_class = :payload,
                optics = optics_prism)
