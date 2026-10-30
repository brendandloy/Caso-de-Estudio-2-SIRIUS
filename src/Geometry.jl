module Geometry

export Debris, generate_debris, generateRandomShape

using LinearAlgebra
using Random
using Statistics
using QHull

const Debris = Dict{Symbol,Any}

# =============================================================================
# Shape generation
# =============================================================================

"""
    generateRandomShape(shapeType, xdim, ydim, zdim)

Returns `(normals, areas, tot_area)` for a triangulated body.

All dimensions are expected in km.
The returned areas are therefore in km^2.

Return:
- normals :: Matrix{Float64} (3 x Nfaces)
- areas   :: Vector{Float64} (Nfaces)
- tot_area:: Float64
"""
function generateRandomShape(shapeType::Int, xdim::Real, ydim::Real, zdim::Real)
    verts = Array{Float64}(undef, 0, 3)
    tris  = Array{Int}(undef, 0, 3)

    if shapeType == 1
        verts, tris = cube_vertices_faces(xdim, ydim, zdim)
    elseif shapeType == 2
        verts, tris = sphere_mesh(xdim, 20, 20)
    elseif shapeType == 3
        verts, tris = cube_vertices_faces(xdim, ydim, zdim)
    elseif shapeType == 4
        verts, tris = cylinder_mesh(xdim, ydim; Nt = 20, Nz = 1)
    elseif shapeType == 5
        verts, tris = hollow_cone(xdim, ydim, false)
    elseif shapeType == 6
        verts, tris = hollow_cone(xdim, ydim, true)
    elseif shapeType == 7
        verts, tris = random_fragment(xdim)
    else
        error("shapeType must be in 1:7.")
    end

    normals, areas, tot_area = calculateNormals(verts, tris)
    return normals, areas, tot_area
end

function cube_vertices_faces(xd, yd, zd)
    x = xd / 2
    y = yd / 2
    z = zd / 2

    V = [
        -x -y -z;
         x -y -z;
         x  y -z;
        -x  y -z;
        -x -y  z;
         x -y  z;
         x  y  z;
        -x  y  z
    ]

    F = [
        1 2 3; 1 3 4;
        5 6 7; 5 7 8;
        1 2 6; 1 6 5;
        2 3 7; 2 7 6;
        3 4 8; 3 8 7;
        4 1 5; 4 5 8
    ]

    return Matrix{Float64}(V), Matrix{Int}(F)
end

"""
    sphere_mesh(radius, nθ, nφ) -> (V, F)

Malla triangular de una esfera centrada en el origen (esfera UV): `nφ - 2` anillos
de latitud con `nθ` vértices cada uno, más un vértice en cada polo unido a su anillo
más cercano con un abanico de triángulos.

Por qué así: si los polos se representan como anillos (como antes), sus `nθ` puntos
coinciden y los triángulos que los tocan tienen área nula y normal NaN.

# Argumentos
- `radius::Real`: radio (en las unidades de la malla; en SIRIUS, km).
- `nθ::Int`: divisiones en longitud (≥ 3).
- `nφ::Int`: niveles de latitud contando los dos polos (≥ 3).

# Devuelve
- `V::Matrix{Float64}` (Nv × 3) con los vértices y `F::Matrix{Int}` (Nf × 3) con los
  triángulos. La orientación de cada normal la fija `calculateNormals`.

# Excepciones
- `ErrorException`: si `nθ < 3` o `nφ < 3`.

# Ejemplo
```julia
V, F = sphere_mesh(1.0, 20, 20)      # 18 anillos + 2 polos: 362 vértices, 720 triángulos
```
"""
function sphere_mesh(radius::Real, nθ::Int, nφ::Int)
    (nθ >= 3 && nφ >= 3) || error("sphere_mesh: nθ and nφ must be >= 3.")

    θs = range(0, 2π, length = nθ + 1)[1:end-1]            # longitudes (sin repetir 2π)
    φs = range(-π / 2, π / 2, length = nφ)[2:end-1]        # latitudes de los anillos (sin los polos)
    nring = length(θs)
    nrings = length(φs)

    # Vértices de los anillos y, al final, polo sur y polo norte
    V = zeros(Float64, nrings * nring + 2, 3)
    idx = (i, j) -> (i - 1) * nring + j                    # vértice j del anillo i
    for (i, φ) in enumerate(φs), (j, θ) in enumerate(θs)
        V[idx(i, j), :] .= (radius * cos(φ) * cos(θ), radius * cos(φ) * sin(θ), radius * sin(φ))
    end
    south = nrings * nring + 1
    north = south + 1
    V[south, :] .= (0.0, 0.0, -radius)
    V[north, :] .= (0.0, 0.0,  radius)

    F = Int[]

    # Bandas entre anillos consecutivos: dos triángulos por cuadrilátero
    for i in 1:(nrings - 1), j in 1:nring
        j2 = (j % nring) + 1
        append!(F, (idx(i, j), idx(i, j2), idx(i + 1, j), idx(i, j2), idx(i + 1, j2), idx(i + 1, j)))
    end

    # Casquetes: abanicos de triángulos desde cada polo a su anillo más cercano
    for j in 1:nring
        j2 = (j % nring) + 1
        append!(F, (south, idx(1, j2), idx(1, j)))
        append!(F, (north, idx(nrings, j), idx(nrings, j2)))
    end

    return V, Matrix{Int}(reshape(F, 3, :)')
end

function cylinder_mesh(R::Real, H::Real; Nt::Int = 20, Nz::Int = 1)
    θs = range(0, 2π, length = Nt + 1)[1:end-1]
    zs = range(-H / 2, H / 2, length = Nz + 2)

    verts = Vector{Tuple{Float64,Float64,Float64}}()
    for z in zs
        for θ in θs
            push!(verts, (R * cos(θ), R * sin(θ), z))
        end
    end

    push!(verts, (0.0, 0.0, -H / 2))
    push!(verts, (0.0, 0.0,  H / 2))

    V = reduce(vcat, (collect(p)' for p in verts))

    nring = length(θs)
    nrings = length(zs)
    idx = (ring, j) -> (ring - 1) * nring + j
    bottom_center = size(V, 1) - 1
    top_center = size(V, 1)

    F = Int[]

    for r in 1:(nrings - 1)
        for j in 1:nring
            j2 = (j % nring) + 1
            v1 = idx(r, j)
            v2 = idx(r, j2)
            v3 = idx(r + 1, j)
            v4 = idx(r + 1, j2)
            append!(F, (v1, v2, v3, v2, v4, v3))
        end
    end

    for j in 1:nring
        j2 = (j % nring) + 1
        append!(F, (bottom_center, idx(1, j2), idx(1, j)))
    end

    for j in 1:nring
        j2 = (j % nring) + 1
        append!(F, (top_center, idx(nrings, j), idx(nrings, j2)))
    end

    F = reshape(F, 3, :)'
    return V, F
end

"""
    hollow_cone(base_radius, height, is_truncated) -> (V, F)

Malla de la superficie lateral de un cono (`is_truncated = false`) o de un tronco de
cono (`is_truncated = true`), con 50 divisiones en el contorno. No incluye bases:
son superficies abiertas (decisión confirmada con el tutor).

Solo se mallan las caras exteriores. Antes había además caras "interiores" con los
mismos vértices; `calculateNormals` las orientaba hacia fuera, así que duplicaban la
superficie (área ×2, ≈ 0.75 mag más brillante). Tampoco tendría sentido orientarlas
hacia dentro: el modelo de Phong no calcula ocultaciones entre caras.

# Argumentos
- `base_radius::Real`: radio de la base mayor.
- `height::Real`: altura.
- `is_truncated::Bool`: si es un tronco; el radio superior es `base_radius / k`, con
  `k` entero aleatorio entre 2 y 10.

# Devuelve
- `V::Matrix{Float64}` (Nv × 3) y `F::Matrix{Int}` (Nf × 3).

# Excepciones
- Ninguna.

# Ejemplo
```julia
V, F = hollow_cone(1.0, 2.0, false)      # 51 vértices (contorno + vértice), 50 triángulos
```
"""
function hollow_cone(base_radius::Real, height::Real, is_truncated::Bool)
    n = 50                                                  # divisiones del contorno
    θ = range(0, 2π, length = n + 1)[1:end-1]

    if is_truncated
        expansion_ratio = rand(2:10)                        # radio mayor / radio menor
        top_radius = base_radius / expansion_ratio

        # Contorno inferior (1..n) y superior (n+1..2n)
        X = vcat(base_radius .* cos.(θ), top_radius .* cos.(θ))
        Y = vcat(base_radius .* sin.(θ), top_radius .* sin.(θ))
        Z = vcat(zeros(n), fill(height, n))

        # Superficie lateral: dos triángulos por segmento del contorno
        F = Int[]
        for i in 1:n
            i2 = (i % n) + 1
            append!(F, (i, i2, n + i2, i, n + i2, n + i))
        end
    else
        # Contorno de la base (1..n) y vértice del cono (n+1)
        X = vcat(base_radius .* cos.(θ), 0.0)
        Y = vcat(base_radius .* sin.(θ), 0.0)
        Z = vcat(zeros(n), height)

        # Superficie lateral: un triángulo por segmento del contorno
        apex = n + 1
        F = Int[]
        for i in 1:n
            i2 = (i % n) + 1
            append!(F, (i, i2, apex))
        end
    end

    return hcat(X, Y, Z), Matrix{Int}(reshape(F, 3, :)')
end

"""
    calculateNormals(V, F) -> (normals, areas, tot_area)

Normales unitarias y áreas de los triángulos de una malla. Cada normal se orienta
hacia fuera respecto al centroide de los vértices (válido para cuerpos convexos o
con forma de estrella respecto a ese punto).

Los triángulos degenerados (vértices alineados o repetidos: área nula y normal
indefinida) se descartan con un aviso, para que no aparezcan normales NaN. Es útil
también para mallas importadas, que a menudo los traen.

# Argumentos
- `V::AbstractMatrix{<:Real}`: vértices (Nv × 3).
- `F::AbstractMatrix{<:Integer}`: triángulos (Nf × 3), índices de `V` empezando en 1.

# Devuelve
- `normals::Matrix{Float64}` (3 × N), `areas::Vector{Float64}` (N) y `tot_area::Float64`,
  con N el número de triángulos no degenerados.

# Excepciones
- `ErrorException`: si todos los triángulos son degenerados.

# Ejemplo
```julia
V, F = cube_vertices_faces(1.0, 1.0, 1.0)
normals, areas, A = calculateNormals(V, F)      # 12 triángulos, A = 6.0
```
"""
function calculateNormals(V::AbstractMatrix{<:Real}, F::AbstractMatrix{<:Integer})
    M = size(F, 1)
    normals = zeros(Float64, 3, M)
    areas = zeros(Float64, M)
    valid = trues(M)                                       # false para triángulos degenerados
    center = vec(mean(V, dims = 1))

    for i in 1:M
        v1, v2, v3 = F[i, 1], F[i, 2], F[i, 3]
        p1 = @view V[v1, :]
        p2 = @view V[v2, :]
        p3 = @view V[v3, :]

        e1 = p2 .- p1
        e2 = p3 .- p1

        n = [
            e1[2] * e2[3] - e1[3] * e2[2],
            e1[3] * e2[1] - e1[1] * e2[3],
            e1[1] * e2[2] - e1[2] * e2[1]
        ]

        # |e1 × e2| despreciable frente a |e1|·|e2|: vértices alineados o repetidos
        if norm(n) <= eps(Float64) * norm(e1) * norm(e2) || iszero(norm(n))
            valid[i] = false
            continue
        end

        Atri = 0.5 * norm(n)
        fc = (p1 .+ p2 .+ p3) ./ 3

        # La normal se gira si apunta hacia el centroide
        if dot(n, fc .- center) < 0
            n .= -n
        end

        normals[:, i] .= n ./ norm(n)
        areas[i] = Atri
    end

    ndeg = count(!, valid)
    ndeg == M && error("calculateNormals: all $M triangles are degenerate.")
    ndeg > 0 && @warn "calculateNormals: $ndeg degenerate triangle(s) (zero area) discarded."

    return normals[:, valid], areas[valid], sum(areas[valid])
end

function generate_points_in_sphere(num_points::Int)
    P = zeros(Float64, num_points, 3)
    count = 0

    while count < num_points
        cand = 2.0 .* rand(3) .- 1.0
        if sum(cand .^ 2) <= 1.0
            count += 1
            P[count, :] .= cand
        end
    end

    return P
end

"""
    random_fragment(char_length) -> (V, F)

Fragmento irregular: envolvente convexa de 10 a 30 puntos aleatorios dentro de una
esfera, escalada para que su longitud característica sea `char_length`.

Longitud característica (definición del NASA Standard Breakup Model):
`L_c = (X + Y + Z) / 3`, con X la dimensión máxima, Y la máxima perpendicular a X y
Z la perpendicular a ambas. Aquí X, Y, Z se miden en los ejes principales de la nube
de puntos (aproximación de esas tres direcciones).

Antes, `char_length` se usaba como si fuera un volumen y un fragmento de 0.1 m salía
con un radio de ≈ 29 m; además la llamada a QHull no existía (`qhull`), así que esta
forma daba error.

# Argumentos
- `char_length::Real`: longitud característica (en las unidades de la malla; en SIRIUS, km).

# Devuelve
- `V::Matrix{Float64}` (Nv × 3), centrado en el origen, y `F::Matrix{Int}` (Nf × 3) con
  los triángulos de la envolvente convexa (los puntos interiores no se usan en `F`).

# Excepciones
- `ErrorException`: si `char_length <= 0`.
- Errores de QHull si los puntos fueran coplanares (no ocurre con puntos aleatorios en 3D).

# Ejemplo
```julia
V, F = random_fragment(1e-4)             # fragmento de L_c = 0.1 m (1e-4 km)
```
"""
function random_fragment(char_length::Real)
    char_length > 0 || error("random_fragment: char_length must be > 0.")

    num_points = rand(10:30)
    pts = generate_points_in_sphere(num_points)            # en la esfera unidad
    pts .-= mean(pts, dims = 1)                            # centrado en su centroide

    # Dimensiones en los ejes principales de la nube de puntos (autovectores de la covarianza)
    pax = eigen(Symmetric(pts' * pts)).vectors
    proj = pts * pax
    extents = vec(maximum(proj, dims = 1) .- minimum(proj, dims = 1))

    # Escala para que (X + Y + Z)/3 = char_length
    pts .*= char_length / mean(extents)

    # Envolvente convexa (QHull vía scipy): simplices son los triángulos, con índices desde 1
    hull = chull(pts)
    return Matrix{Float64}(pts), Matrix{Int}(hull.simplices)
end

# =============================================================================
# Debris generation
# =============================================================================

function _rand_trunc_normal(μ, σ; upper = 1.0, lower = 0.0)
    x = μ + σ * randn()
    while (x < lower) || (x >= upper)
        x = μ + σ * randn()
    end
    return x
end

# Normal truncada a valores positivos: dimensiones aleatorias sin valores negativos
_pos_normal(μ, σ) = _rand_trunc_normal(μ, σ; lower = 0.0, upper = Inf)

function _make_optical_coeffs_body(μ_abs, σ_abs, μ_spec_aux, σ_spec_aux)
    Cabs = _rand_trunc_normal(μ_abs, σ_abs)
    Cesp_aux = _rand_trunc_normal(μ_spec_aux, σ_spec_aux)
    Cesp = Cesp_aux * (1 - Cabs)
    Cdif = max(1 - Cabs - Cesp, 0.0)
    return Cabs, Cesp, Cdif
end

function _dims_prism(sat_size::Int, sub_cat::Int)
    xdim = 0.0
    ydim = 0.0
    zdim = 0.0
    siz = 0
    siz2 = 0.0

    if sat_size == 1
        if sub_cat == 1
            xdim = 0.1; ydim = xdim; zdim = xdim; siz2 = 1.1
        elseif sub_cat == 2
            xdim = 0.1; ydim = xdim; zdim = 3 * xdim; siz2 = 1.2
        elseif sub_cat == 3
            xdim = 0.1; ydim = 2 * xdim; zdim = 3 * xdim; siz2 = 1.3
        else
            # xdim = 0.2; ydim = xdim; zdim = 1.5 * xdim; siz2 = 1.4
            xdim = 4.1; ydim = 4.1; zdim = 2.7; siz2 = 1.4
        end
        siz = 1
    elseif sat_size == 2
        if sub_cat == 1
            xdim = _pos_normal(0.5, 0.1)
            ydim = _pos_normal(0.5, 0.1)
            zdim = _pos_normal(0.5, 0.1)
            siz2 = 2.1
        elseif sub_cat == 2
            xdim = _pos_normal(0.7, 0.1)
            ydim = _pos_normal(0.7, 0.1)
            zdim = _pos_normal(0.7, 0.1)
            siz2 = 2.2
        elseif sub_cat == 3
            xdim = _pos_normal(1.0, 0.1)
            ydim = _pos_normal(1.0, 0.1)
            zdim = _pos_normal(1.0, 0.1)
            siz2 = 2.3
        else
            xdim = 1.1; ydim = 1.5; zdim = 1.2; siz2 = 2.4
        end
        siz = 2
    else
        if sub_cat == 1
            xdim = _pos_normal(1.25, 0.15)
            ydim = _pos_normal(1.25, 0.15)
            zdim = _pos_normal(1.25, 0.15)
            siz2 = 3.1
        elseif sub_cat == 2
            xdim = _pos_normal(1.5, 0.15)
            ydim = _pos_normal(1.5, 0.15)
            zdim = _pos_normal(1.5, 0.15)
            siz2 = 3.2
        elseif sub_cat == 3
            xdim = _pos_normal(1.8, 0.15)
            ydim = _pos_normal(1.8, 0.15)
            zdim = _pos_normal(1.8, 0.15)
            siz2 = 3.3
        else
            xdim = 4.5; ydim = 4.5; zdim = 5.0; siz2 = 3.4
        end
        siz = 3
    end

    return xdim, ydim, zdim, siz, siz2
end

function _dims_sphere(sat_size::Int)
    if sat_size == 1
        return 0.1, 1, 1.1
    elseif sat_size == 2
        return 0.3, 2, 2.1
    else
        return 0.5, 3, 3.1
    end
end

function _dims_plate(sat_size::Int, sub_cat::Int)
    xdim = 0.0
    ydim = 0.0
    zdim = 0.0
    siz = 0
    siz2 = 0.0

    if sat_size == 1
        if sub_cat == 1
            xdim = 0.1; ydim = xdim; zdim = 0.001; siz2 = 1.1
        elseif sub_cat == 2
            xdim = 0.2; ydim = 2 * xdim; zdim = 0.002; siz2 = 1.2
        elseif sub_cat == 3
            xdim = 0.3; ydim = 3 * xdim; zdim = 0.003; siz2 = 1.3
        else
            xdim = 0.1 * rand(1:5); ydim = 4 * xdim; zdim = 0.04 * xdim; siz2 = 1.4
        end
        siz = 1
    elseif sat_size == 2
        if sub_cat == 1
            xdim = 0.5; ydim = 2 * xdim; zdim = 0.005; siz2 = 2.1
        elseif sub_cat == 2
            xdim = 0.6; ydim = 2 * xdim; zdim = 0.006; siz2 = 2.2
        elseif sub_cat == 3
            xdim = 0.7; ydim = 2 * xdim; zdim = 0.007; siz2 = 2.3
        else
            xdim = 0.1 * rand(5:10); ydim = 2 * xdim; zdim = 0.08 * xdim; siz2 = 2.4
        end
        siz = 2
    else
        if sub_cat == 1
            xdim = 1.0; ydim = 2 * xdim; zdim = 0.01; siz2 = 3.1
        elseif sub_cat == 2
            xdim = 1.5; ydim = 2 * xdim; zdim = 0.01; siz2 = 3.2
        elseif sub_cat == 3
            xdim = 2.0; ydim = 2 * xdim; zdim = 0.02; siz2 = 3.3
        else
            xdim = 2.5; ydim = 2 * xdim; zdim = 0.05; siz2 = 3.4
        end
        siz = 3
    end

    return xdim, ydim, zdim, siz, siz2
end

function _dims_cylinder(sat_size::Int, sub_cat::Int)
    radius = 0.0
    height = 0.0
    siz = 0
    siz2 = 0.0

    if sat_size == 1
        if sub_cat == 1
            radius = 0.1; height = 0.3; siz2 = 1.1
        elseif sub_cat == 2
            radius = 0.15; height = 0.5; siz2 = 1.2
        elseif sub_cat == 3
            radius = 0.2; height = 0.6; siz2 = 1.3
        else
            radius = 0.3; height = 0.8; siz2 = 1.4
        end
        siz = 1
    elseif sat_size == 2
        if sub_cat == 1
            radius = 0.3; height = 1.0; siz2 = 2.1
        elseif sub_cat == 2
            radius = 0.4; height = 1.5; siz2 = 2.2
        elseif sub_cat == 3
            radius = 0.5; height = 2.0; siz2 = 2.3
        else
            radius = 0.6; height = 2.5; siz2 = 2.4
        end
        siz = 2
    else
        if sub_cat == 1
            radius = 0.8; height = 3.0; siz2 = 3.1
        elseif sub_cat == 2
            radius = 1.0; height = 3.5; siz2 = 3.2
        elseif sub_cat == 3
            radius = 1.2; height = 4.0; siz2 = 3.3
        else
            radius = 1.5; height = 5.0; siz2 = 3.4
        end
        siz = 3
    end

    return radius, height, siz, siz2
end

function _dims_cone(sat_size::Int, sub_cat::Int)
    ratio = 0.1 * rand(1:5)
    xdim = 0.0
    ydim = 0.0
    zdim = 0.0
    siz = 0
    siz2 = 0.0

    if sat_size == 1
        if sub_cat == 1
            xdim = _pos_normal(0.1, 0.05)
        elseif sub_cat == 2
            xdim = _pos_normal(0.2, 0.05)
        elseif sub_cat == 3
            xdim = _pos_normal(0.3, 0.05)
        else
            xdim = 0.1 * rand(1:5)
        end
        siz = 1
        siz2 = 1.0 + 0.1 * sub_cat
    elseif sat_size == 2
        if sub_cat == 1
            xdim = _pos_normal(0.6, 0.05)
        elseif sub_cat == 2
            xdim = _pos_normal(0.7, 0.05)
        elseif sub_cat == 3
            xdim = _pos_normal(0.8, 0.05)
        else
            xdim = _pos_normal(0.9, 0.05)
        end
        siz = 2
        siz2 = 2.0 + 0.1 * sub_cat
    else
        if sub_cat == 1
            xdim = _pos_normal(1.0, 0.05)
        elseif sub_cat == 2
            xdim = _pos_normal(1.5, 0.05)
        elseif sub_cat == 3
            xdim = _pos_normal(2.0, 0.05)
        else
            xdim = _pos_normal(2.5, 0.05)
        end
        siz = 3
        siz2 = 3.0 + 0.1 * sub_cat
    end

    ydim = xdim / ratio
    zdim = xdim / 2
    return xdim, ydim, zdim, siz, siz2
end

function _dims_trunk(sat_size::Int, sub_cat::Int)
    ratio = rand(1:5)
    xdim = 0.0
    ydim = 0.0
    zdim = 0.0
    siz = 0
    siz2 = 0.0

    if sat_size == 1
        if sub_cat == 1
            xdim = _pos_normal(0.1, 0.05)
        elseif sub_cat == 2
            xdim = _pos_normal(0.2, 0.05)
        elseif sub_cat == 3
            xdim = _pos_normal(0.3, 0.05)
        else
            xdim = 0.1 * rand(1:5)
        end
        siz = 1
        siz2 = 1.0 + 0.1 * sub_cat
    elseif sat_size == 2
        if sub_cat == 1
            xdim = _pos_normal(0.6, 0.05)
        elseif sub_cat == 2
            xdim = _pos_normal(0.7, 0.05)
        elseif sub_cat == 3
            xdim = _pos_normal(0.8, 0.05)
        else
            xdim = _pos_normal(0.9, 0.05)
        end
        siz = 2
        siz2 = 2.0 + 0.1 * sub_cat
    else
        if sub_cat == 1
            xdim = _pos_normal(1.0, 0.05)
        elseif sub_cat == 2
            xdim = _pos_normal(1.5, 0.05)
        elseif sub_cat == 3
            xdim = _pos_normal(2.0, 0.05)
        else
            xdim = _pos_normal(2.5, 0.05)
        end
        siz = 3
        siz2 = 3.0 + 0.1 * sub_cat
    end

    ydim = xdim * ratio
    zdim = xdim / 2
    return xdim, ydim, zdim, siz, siz2
end

function _dims_fragment(sat_size::Int)
    if sat_size == 1
        return 0.1, 1, 1.1
    elseif sat_size == 2
        return 0.3, 2, 2.1
    else
        return 0.5, 3, 3.1
    end
end


function generate_debris(sat_size::Int, sub_cat::Int;
                         deploy_panels::Bool = false,
                         shape_type::Int = 1)
    debris = Debris()

    debris[:sat_size] = sat_size
    debris[:sub_cat] = sub_cat
    debris[:deploy_panels] = deploy_panels
    debris[:shape_type] = shape_type
    debris[:alpha] = max(0.1, 1.2 + 0.1 * randn())

    Cabs_panel, Cesp_panel, Cdif_panel = _make_optical_coeffs_body(0.75, 0.1, 0.8, 0.05)

    if deploy_panels
        Cabs_body, Cesp_body, Cdif_body = _make_optical_coeffs_body(0.62, 0.12, 0.7, 0.1)

        xdim_m, ydim_m, zdim_m, siz, siz2 = _dims_prism(sat_size, sub_cat)

        # xdim = 1e-3 *(xdim_m)
        # ydim = 1e-3 *(ydim_m)
        # zdim = 1e-3 *(zdim_m)

        xdim = 1e-3 *(4.1)
        ydim = 1e-3 *(4.1)
        zdim = 1e-3 *(12.8)

        normals, areas, total_area = generateRandomShape(1, xdim, ydim, zdim)

        debris[:shape_name] = :prism_with_panels
        debris[:siz] = siz
        debris[:siz2] = siz2
        debris[:normals] = normals
        debris[:areas] = areas
        debris[:total_area] = total_area
        debris[:dims] = (xdim = xdim, ydim = ydim, zdim = zdim)

        A_panels = 2 * xdim * ydim
        A_body = 2 * (xdim * zdim + ydim * zdim)
        A_geom = A_panels + A_body

        debris[:C_abs_body] = Cabs_body
        debris[:C_spec_body] = Cesp_body
        debris[:C_diff_body] = Cdif_body
        debris[:C_abs_panel] = Cabs_panel
        debris[:C_spec_panel] = Cesp_panel
        debris[:C_diff_panel] = Cdif_panel

        debris[:C_spec] = (Cesp_panel * A_panels + Cesp_body * A_body) / A_geom
        debris[:C_diff] = (Cdif_panel * A_panels + Cdif_body * A_body) / A_geom
        debris[:C_abs] = 1 - debris[:C_spec] - debris[:C_diff]

        return debris
    end

    Cabs_body, Cesp_body, Cdif_body = _make_optical_coeffs_body(0.62, 0.12, 0.7, 0.1)

    if shape_type == 1
        Cabs_top, Cesp_top, Cdif_top = _make_optical_coeffs_body(0.55, 0.2, 0.7, 0.1)

        xdim_m, ydim_m, zdim_m, siz, siz2 = _dims_prism(sat_size, sub_cat)

        xdim = 1e-3 *(xdim_m)
        ydim = 1e-3 *(ydim_m)
        zdim = 1e-3 *(zdim_m)

        normals, areas, total_area = generateRandomShape(1, xdim, ydim, zdim)

        debris[:shape_name] = :prism_no_panels
        debris[:siz] = siz
        debris[:siz2] = siz2
        debris[:normals] = normals
        debris[:areas] = areas
        debris[:total_area] = total_area
        debris[:dims] = (xdim = xdim, ydim = ydim, zdim = zdim)

        A_top = 2 * xdim * ydim
        A_body = 2 * (xdim * zdim + ydim * zdim)
        A_geom = A_top + A_body

        debris[:C_abs_body] = Cabs_body
        debris[:C_spec_body] = Cesp_body
        debris[:C_diff_body] = Cdif_body
        debris[:C_abs_top] = Cabs_top
        debris[:C_spec_top] = Cesp_top
        debris[:C_diff_top] = Cdif_top

        debris[:C_spec] = (Cesp_top * A_top + Cesp_body * A_body) / A_geom
        debris[:C_diff] = (Cdif_top * A_top + Cdif_body * A_body) / A_geom
        debris[:C_abs] = 1 - debris[:C_spec] - debris[:C_diff]

    elseif shape_type == 2
        radius_m, siz, siz2 = _dims_sphere(sat_size)
        radius = 1e-3 *(radius_m)

        normals, areas, total_area = generateRandomShape(2, radius, 0.0, 0.0)

        debris[:shape_name] = :sphere
        debris[:siz] = siz
        debris[:siz2] = siz2
        debris[:normals] = normals
        debris[:areas] = areas
        debris[:total_area] = total_area
        debris[:dims] = (radius = radius,)
        debris[:C_abs] = Cabs_body
        debris[:C_spec] = Cesp_body
        debris[:C_diff] = Cdif_body

    elseif shape_type == 3
        xdim_m, ydim_m, zdim_m, siz, siz2 = _dims_plate(sat_size, sub_cat)

        xdim = 1e-3 *(xdim_m)
        ydim = 1e-3 *(ydim_m)
        zdim = 1e-3 *(zdim_m)

        normals, areas, total_area = generateRandomShape(3, xdim, ydim, zdim)

        debris[:shape_name] = :plate
        debris[:siz] = siz
        debris[:siz2] = siz2
        debris[:normals] = normals
        debris[:areas] = areas
        debris[:total_area] = total_area
        debris[:dims] = (xdim = xdim, ydim = ydim, zdim = zdim)
        debris[:C_abs] = Cabs_body
        debris[:C_spec] = Cesp_panel
        debris[:C_diff] = Cdif_panel

    elseif shape_type == 4
        radius_m, height_m, siz, siz2 = _dims_cylinder(sat_size, sub_cat)

        radius = 1e-3 *(radius_m)
        height = 1e-3 *(height_m)

        normals, areas, total_area = generateRandomShape(4, radius, height, 0.0)

        debris[:shape_name] = :cylinder
        debris[:siz] = siz
        debris[:siz2] = siz2
        debris[:normals] = normals
        debris[:areas] = areas
        debris[:total_area] = total_area
        debris[:dims] = (radius = radius, height = height)
        debris[:C_abs] = Cabs_body
        debris[:C_spec] = Cesp_body
        debris[:C_diff] = Cdif_body

    elseif shape_type == 5
        xdim_m, ydim_m, zdim_m, siz, siz2 = _dims_cone(sat_size, sub_cat)

        xdim = 1e-3 *(xdim_m)
        ydim = 1e-3 *(ydim_m)
        zdim = 1e-3 *(zdim_m)

        normals, areas, total_area = generateRandomShape(5, xdim, ydim, zdim)

        debris[:shape_name] = :cone
        debris[:siz] = siz
        debris[:siz2] = siz2
        debris[:normals] = normals
        debris[:areas] = areas
        debris[:total_area] = total_area
        debris[:dims] = (xdim = xdim, ydim = ydim, zdim = zdim)
        debris[:C_abs] = Cabs_body
        debris[:C_spec] = Cesp_body
        debris[:C_diff] = Cdif_body

    elseif shape_type == 6
        xdim_m, ydim_m, zdim_m, siz, siz2 = _dims_trunk(sat_size, sub_cat)

        xdim = 1e-3 *(xdim_m)
        ydim = 1e-3 *(ydim_m)
        zdim = 1e-3 *(zdim_m)

        normals, areas, total_area = generateRandomShape(6, xdim, ydim, zdim)

        debris[:shape_name] = :trunk
        debris[:siz] = siz
        debris[:siz2] = siz2
        debris[:normals] = normals
        debris[:areas] = areas
        debris[:total_area] = total_area
        debris[:dims] = (xdim = xdim, ydim = ydim, zdim = zdim)
        debris[:C_abs] = Cabs_body
        debris[:C_spec] = Cesp_body
        debris[:C_diff] = Cdif_body

    elseif shape_type == 7
        xdim_m, siz, siz2 = _dims_fragment(sat_size)
        xdim = 1e-3 *(xdim_m)

        normals, areas, total_area = generateRandomShape(7, xdim, 0.0, 0.0)

        debris[:shape_name] = :fragment
        debris[:siz] = siz
        debris[:siz2] = siz2
        debris[:normals] = normals
        debris[:areas] = areas
        debris[:total_area] = total_area
        debris[:dims] = (char_length = xdim,)
        debris[:C_abs] = Cabs_body
        debris[:C_spec] = Cesp_body
        debris[:C_diff] = Cdif_body

    else
        error("shape_type must be in 1:7.")
    end

    return debris
end

end # module Geometry