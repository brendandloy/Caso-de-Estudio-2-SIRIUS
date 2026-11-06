# =============================================================================
# Meshes.jl
#
# Mallas triangulares de las formas de SIRIUS y cálculo de normales y áreas.
# Se incluye desde Geometry.jl, por lo que comparte su espacio de nombres.
#
# Unidades: las de las dimensiones de entrada (en SIRIUS, km; las áreas salen en km²).
# Ejes cuerpo: los de cada función (ver la tabla de ejes en Geometry.jl).
#
# Funciones:
#   cube_vertices_faces       : prisma rectangular (prisma y placa)
#   sphere_mesh               : esfera UV con un vértice en cada polo
#   cylinder_mesh             : cilindro cerrado (lateral + dos tapas)
#   hollow_cone               : superficie lateral de un cono o de un tronco de cono
#   calculateNormals          : normales unitarias hacia fuera y áreas de los triángulos
#   generate_points_in_sphere : puntos uniformes dentro de la esfera unidad
#   convex_hull               : envolvente convexa de una nube de puntos pequeña
#   principal_axes            : rota una nube de puntos a sus ejes principales
#   random_fragment           : fragmento irregular con longitud característica dada
# =============================================================================

"""
    cube_vertices_faces(xd, yd, zd) -> (V, F)

Malla de un prisma rectangular centrado en el origen, con las aristas paralelas a
los ejes x, y, z: 8 vértices y 12 triángulos (dos por cara).

# Argumentos
- `xd`, `yd`, `zd`: longitud de las aristas según x, y, z.

# Devuelve
- `V::Matrix{Float64}` (8 × 3) y `F::Matrix{Int}` (12 × 3).

# Excepciones
- Ninguna.

# Ejemplo
```julia
V, F = cube_vertices_faces(1.0, 2.0, 3.0)
```
"""
function cube_vertices_faces(xd, yd, zd)
    x = xd / 2
    y = yd / 2
    z = zd / 2

    # Vértices: cara inferior (z = -z) y superior (z = +z)
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

    # Dos triángulos por cara: -z, +z, -y, +x, +y, -x
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

Por qué así: si los polos se representan como anillos, sus `nθ` puntos coinciden y
los triángulos que los tocan tienen área nula y normal NaN.

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

"""
    cylinder_mesh(R, H; Nt=20, Nz=1) -> (V, F)

Malla de un cilindro cerrado centrado en el origen, con eje según z (de `-H/2` a
`+H/2`): superficie lateral en `Nz + 1` bandas y dos tapas en abanico.

# Argumentos
- `R::Real`: radio.
- `H::Real`: altura.
- `Nt::Int`: divisiones del contorno (por defecto 20).
- `Nz::Int`: anillos intermedios de la superficie lateral (por defecto 1).

# Devuelve
- `V::Matrix{Float64}` (Nv × 3) y `F` (Nf × 3) con los triángulos.

# Excepciones
- Ninguna.

# Ejemplo
```julia
V, F = cylinder_mesh(1.0, 3.0)        # 62 vértices, 120 triángulos
```
"""
function cylinder_mesh(R::Real, H::Real; Nt::Int = 20, Nz::Int = 1)
    θs = range(0, 2π, length = Nt + 1)[1:end-1]            # ángulos del contorno
    zs = range(-H / 2, H / 2, length = Nz + 2)             # alturas de los anillos

    # Anillos del contorno, de abajo arriba
    verts = Vector{Tuple{Float64,Float64,Float64}}()
    for z in zs
        for θ in θs
            push!(verts, (R * cos(θ), R * sin(θ), z))
        end
    end

    # Centros de las tapas inferior y superior
    push!(verts, (0.0, 0.0, -H / 2))
    push!(verts, (0.0, 0.0,  H / 2))

    V = reduce(vcat, (collect(p)' for p in verts))

    nring = length(θs)
    nrings = length(zs)
    idx = (ring, j) -> (ring - 1) * nring + j               # vértice j del anillo ring
    bottom_center = size(V, 1) - 1
    top_center = size(V, 1)

    F = Int[]

    # Superficie lateral: dos triángulos por cuadrilátero
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

    # Tapa inferior
    for j in 1:nring
        j2 = (j % nring) + 1
        append!(F, (bottom_center, idx(1, j2), idx(1, j)))
    end

    # Tapa superior
    for j in 1:nring
        j2 = (j % nring) + 1
        append!(F, (top_center, idx(nrings, j), idx(nrings, j2)))
    end

    F = reshape(F, 3, :)'
    return V, F
end

"""
    hollow_cone(base_radius, height; top_radius=0.0) -> (V, F)

Malla de la superficie lateral de un cono (`top_radius = 0`) o de un tronco de cono
(`0 < top_radius < base_radius`), con 50 divisiones en el contorno. La base mayor
está en `z = 0` y la menor (o el vértice) en `z = height`. No incluye bases: son
superficies abiertas (decisión confirmada con el tutor).

Solo se mallan las caras exteriores: unas caras "interiores" con los mismos vértices
duplicarían la superficie, porque `calculateNormals` las orienta hacia fuera, y el
modelo de Phong no calcula ocultaciones entre caras.

Antes el radio superior del tronco se sorteaba aquí dentro (base / k, k entre 2 y 10);
ahora es un parámetro explícito.

# Argumentos
- `base_radius::Real`: radio de la base mayor.
- `height::Real`: altura.
- `top_radius::Real`: radio de la base menor (0 para un cono).

# Devuelve
- `V::Matrix{Float64}` (Nv × 3) y `F::Matrix{Int}` (Nf × 3).

# Excepciones
- `ErrorException`: si `top_radius` es negativo o no es menor que `base_radius`.

# Ejemplo
```julia
V, F = hollow_cone(1.0, 2.0)                      # cono: 51 vértices, 50 triángulos
V, F = hollow_cone(1.0, 2.0; top_radius = 0.5)    # tronco: 100 vértices, 100 triángulos
```
"""
function hollow_cone(base_radius::Real, height::Real; top_radius::Real = 0.0)
    (0 <= top_radius < base_radius) ||
        error("hollow_cone: top_radius must be in [0, base_radius) (got $top_radius).")

    n = 50                                                  # divisiones del contorno
    θ = range(0, 2π, length = n + 1)[1:end-1]

    if top_radius > 0
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
indefinida) se descartan con un aviso, para que no aparezcan normales NaN.

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
    is_valid = trues(M)                                       # false para triángulos degenerados
    center = vec(mean(V, dims = 1))

    for i in 1:M
        v1, v2, v3 = F[i, 1], F[i, 2], F[i, 3]
        p1 = @view V[v1, :]
        p2 = @view V[v2, :]
        p3 = @view V[v3, :]

        e1 = p2 .- p1
        e2 = p3 .- p1

        # Producto vectorial e1 × e2: normal sin normalizar, de módulo 2·área
        n = [
            e1[2] * e2[3] - e1[3] * e2[2],
            e1[3] * e2[1] - e1[1] * e2[3],
            e1[1] * e2[2] - e1[2] * e2[1]
        ]

        # |e1 × e2| despreciable frente a |e1|·|e2|: vértices alineados o repetidos
        if norm(n) <= eps(Float64) * norm(e1) * norm(e2) || iszero(norm(n))
            is_valid[i] = false
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

    ndeg = count(!, is_valid)
    ndeg == M && error("calculateNormals: all $M triangles are degenerate.")
    ndeg > 0 && @warn "calculateNormals: $ndeg degenerate triangle(s) (zero area) discarded."

    return normals[:, is_valid], areas[is_valid], sum(areas[is_valid])
end

"""
    generate_points_in_sphere(rng, num_points) -> Matrix{Float64}

Puntos uniformes dentro de la esfera unidad, por rechazo dentro del cubo [-1, 1]³.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio (con la semilla del objetivo).
- `num_points::Int`: número de puntos.

# Devuelve
- `Matrix{Float64}` (num_points × 3).

# Excepciones
- Ninguna.
"""
function generate_points_in_sphere(rng::AbstractRNG, num_points::Int)
    P = zeros(Float64, num_points, 3)
    count = 0

    while count < num_points
        cand = 2.0 .* rand(rng, 3) .- 1.0                  # punto uniforme en el cubo
        if sum(cand .^ 2) <= 1.0                           # se acepta si cae dentro de la esfera
            count += 1
            P[count, :] .= cand
        end
    end

    return P
end

"""
    convex_hull(P) -> Matrix{Int}

Envolvente convexa de una nube pequeña de puntos en posición general, por fuerza
bruta: un triángulo (i, j, k) es cara de la envolvente si todos los demás puntos
quedan a un mismo lado de su plano. Con N ≤ 30 puntos son unas 4000 ternas.

Por qué existe: sustituye a QHull (que necesitaba scipy a través de PyCall) para que
Geometry no tenga dependencias externas y la pueda cargar también el generador de JSON.
Con puntos aleatorios no hay cuatro coplanares (probabilidad nula), así que cada cara
es un triángulo y no se duplica.

# Argumentos
- `P::AbstractMatrix{<:Real}`: puntos (N × 3), N ≥ 4, no coplanares.

# Devuelve
- `Matrix{Int}` (Nf × 3): triángulos de la envolvente (la orientación la fija `calculateNormals`).

# Excepciones
- `ErrorException`: si hay menos de 4 puntos o no se encuentra ninguna cara (puntos coplanares).

# Ejemplo
```julia
P = [0.0 0 0; 1 0 0; 0 1 0; 0 0 1; 0.1 0.1 0.1]
convex_hull(P)          # 4 caras: el punto interior no aparece
```
"""
function convex_hull(P::AbstractMatrix{<:Real})
    N = size(P, 1)
    N >= 4 || error("convex_hull: at least 4 points are needed (got $N).")
    scale = maximum(abs, P)                                # tamaño de la nube, para la tolerancia

    F = Int[]
    for i in 1:N-2, j in i+1:N-1, k in j+1:N
        nrm = cross(P[j, :] .- P[i, :], P[k, :] .- P[i, :])   # normal del plano (i, j, k)
        norm(nrm) > 0 || continue
        d = (P .- P[i:i, :]) * nrm                         # distancia (con signo) × |nrm| al plano
        tol = 1e-10 * norm(nrm) * scale                    # tolerancia relativa al tamaño
        d[[i, j, k]] .= 0.0
        # Cara de la envolvente: ningún punto a los dos lados del plano
        (all(<=(tol), d) || all(>=(-tol), d)) && append!(F, (i, j, k))
    end

    isempty(F) && error("convex_hull: no faces found (coplanar points?).")
    return Matrix{Int}(reshape(F, 3, :)')
end

"""
    principal_axes(P) -> (Pr, extents)

Expresa una nube de puntos centrada en sus ejes principales (autovectores de la
matriz de covarianza), ordenados de mayor a menor extensión: x es la dimensión
mayor, y la siguiente y z la menor. La rotación es propia (det = +1).

# Argumentos
- `P::AbstractMatrix{<:Real}`: puntos (N × 3), ya centrados en su centroide.

# Devuelve
- `Pr::Matrix{Float64}` (N × 3): puntos en los ejes principales.
- `extents::Vector{Float64}`: extensión (máx - mín) según x, y, z, de mayor a menor.

# Excepciones
- Ninguna.
"""
function principal_axes(P::AbstractMatrix{<:Real})
    vecs = eigen(Symmetric(P' * P)).vectors                # columnas: ejes principales
    proj = P * vecs
    ext = vec(maximum(proj, dims = 1) .- minimum(proj, dims = 1))

    # Columnas ordenadas de mayor a menor extensión; se corrige el signo si sale una reflexión
    perm = sortperm(ext; rev = true)
    R = vecs[:, perm]
    det(R) < 0 && (R[:, 3] .*= -1)

    Pr = P * R
    return Pr, ext[perm]
end

"""
    random_fragment(char_length; rng=Random.default_rng()) -> (V, F)

Fragmento irregular: envolvente convexa de 10 a 30 puntos aleatorios dentro de una
esfera, expresada en sus ejes principales y escalada para que su longitud
característica sea `char_length`.

Longitud característica (definición del NASA Standard Breakup Model):
`L_c = (X + Y + Z) / 3`, con X la dimensión máxima, Y la máxima perpendicular a X y
Z la perpendicular a ambas. Aquí X, Y, Z son las extensiones en los ejes principales
de la nube, que son los ejes cuerpo del fragmento (x = X, y = Y, z = Z).

# Argumentos
- `char_length::Real`: longitud característica (en las unidades de la malla; en SIRIUS, km).
- `rng::AbstractRNG`: generador aleatorio (con la semilla del objetivo, para poder repetirlo).

# Devuelve
- `V::Matrix{Float64}` (Nv × 3), centrado en el origen, y `F::Matrix{Int}` (Nf × 3) con
  los triángulos de la envolvente convexa (los puntos interiores no se usan en `F`).

# Excepciones
- `ErrorException`: si `char_length <= 0`.

# Ejemplo
```julia
V, F = random_fragment(1e-4; rng = MersenneTwister(1))     # fragmento de L_c = 0.1 m (1e-4 km)
```
"""
function random_fragment(char_length::Real; rng::AbstractRNG = Random.default_rng())
    char_length > 0 || error("random_fragment: char_length must be > 0.")

    num_points = rand(rng, 10:30)
    pts = generate_points_in_sphere(rng, num_points)       # en la esfera unidad
    pts .-= mean(pts, dims = 1)                            # centrado en su centroide

    # Ejes cuerpo = ejes principales (x mayor, z menor)
    pts, extents = principal_axes(pts)

    # Escala para que (X + Y + Z)/3 = char_length
    pts .*= char_length / mean(extents)

    return pts, convex_hull(pts)
end
