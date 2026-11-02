# =============================================================================
# plate.jl
#
# Forma "plate": placa rectangular delgada (tapas, paneles sueltos, cofias aplanadas).
#
#   "shape": {"type": "plate", "x_m": 1.0, "y_m": 0.5, "thickness_m": 0.01}
#
# Ejes cuerpo: origen en el centro; x, y en el plano de la placa y z normal a ella
# (thickness_m), que tiene que ser la dimensión menor. Material de panel.
# Sorteado: x ≥ y, x/y en PLATE_ASPECT_RANGE y t/x en PLATE_THICKNESS_RATIO_RANGE.
# Sin medidas, el tamaño sale de los objetos de misión. L_c = (x + y + t)/3.
#
# Funciones: check_plate, unit_plate, lc_plate, surface_plate, mesh_plate, optics_plate
# =============================================================================

"""
    check_plate(d, who)

Comprueba que el espesor es la dimensión menor de la placa (si no, sería un prisma).

# Argumentos
- `d`: medidas (m).
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si `thickness_m` es mayor que `x_m` o que `y_m`.
"""
function check_plate(d, who::AbstractString)
    d["thickness_m"] <= min(d["x_m"], d["y_m"]) ||
        error("$who: 'thickness_m' must be the smallest dimension (use 'prism' otherwise).")
    return nothing
end

"""
    unit_plate(rng) -> Dict{String,Float64}

Medidas de la placa a escala 1: x = 1, y = 1/(x/y) con x/y uniforme en
`PLATE_ASPECT_RANGE` y espesor uniforme en `PLATE_THICKNESS_RATIO_RANGE`.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.

# Devuelve
- `Dict` con `x_m`, `y_m`, `thickness_m`.

# Excepciones
- Ninguna.
"""
function unit_plate(rng::AbstractRNG)
    alo, ahi = PLATE_ASPECT_RANGE
    tlo, thi = PLATE_THICKNESS_RATIO_RANGE
    return Dict{String,Float64}("x_m" => 1.0,
                                "y_m" => 1 / (alo + (ahi - alo) * rand(rng)),
                                "thickness_m" => tlo + (thi - tlo) * rand(rng))
end

# L_c = (x + y + t)/3
lc_plate(d) = (d["x_m"] + d["y_m"] + d["thickness_m"]) / 3

# Superficie 2(xy + xt + yt)
surface_plate(d) = 2 * (d["x_m"] * d["y_m"] + (d["x_m"] + d["y_m"]) * d["thickness_m"])

# Malla de 12 triángulos con el espesor según z (medidas en km)
mesh_plate(d, rng::AbstractRNG) = cube_vertices_faces(d["x_m"], d["y_m"], d["thickness_m"])

# Material de panel en todas las caras
optics_plate(rng::AbstractRNG, d) = optical_coeffs(rng, OPTICS_PANEL)

register_shape!("plate"; dims = ("x_m", "y_m", "thickness_m"), unit = unit_plate, lc = lc_plate,
                surface = surface_plate, mesh = mesh_plate, size_class = :mission_related,
                check = check_plate, optics = optics_plate)
