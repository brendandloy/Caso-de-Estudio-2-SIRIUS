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

function sphere_mesh(radius::Real, nθ::Int, nφ::Int)
    θs = range(0, 2π, length = nθ + 1)[1:end-1]
    φs = range(-π / 2, π / 2, length = nφ)

    points = Vector{NTuple{3,Float64}}()
    for φ in φs, θ in θs
        x = radius * cos(φ) * cos(θ)
        y = radius * cos(φ) * sin(θ)
        z = radius * sin(φ)
        push!(points, (x, y, z))
    end

    V = reduce(vcat, (collect(p)' for p in points))

    F = Int[]
    idx = (i, j) -> (i - 1) * length(θs) + j

    for i in 1:(length(φs) - 1)
        for j in 1:length(θs)
            j2 = (j % length(θs)) + 1
            v1 = idx(i, j)
            v2 = idx(i, j2)
            v3 = idx(i + 1, j)
            v4 = idx(i + 1, j2)
            append!(F, (v1, v2, v3, v2, v4, v3))
        end
    end

    F = reshape(F, 3, :)'
    return V, F
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

function hollow_cone(base_radius::Real, height::Real, is_truncated::Bool)
    n = 50
    θ = range(0, 2π, length = n + 1)[1:end-1]

    if is_truncated
        expansion_ratio = rand(2:10)
        top_radius = base_radius / expansion_ratio

        X = vcat(base_radius .* cos.(θ), top_radius .* cos.(θ))
        Y = vcat(base_radius .* sin.(θ), top_radius .* sin.(θ))
        Z = vcat(zeros(n), fill(height, n))

        Fext = Int[]
        for i in 1:n
            i2 = (i % n) + 1
            append!(Fext, (i, i2, n + i2, i, n + i2, n + i))
        end

        Fint = Int[]
        for i in 1:n
            i2 = (i % n) + 1
            append!(Fint, (n + i2, n + i, i, n + i2, i, i2))
        end

        V = hcat(X, Y, Z)
        F = reshape(vcat(Fext, Fint), 3, :)'
        return V, F
    else
        X = vcat(base_radius .* cos.(θ), 0.0)
        Y = vcat(base_radius .* sin.(θ), 0.0)
        Z = vcat(zeros(n), height)

        apex = n + 1
        Fext = Int[]
        Fint = Int[]

        for i in 1:n
            i2 = (i % n) + 1
            append!(Fext, (i, i2, apex))
            append!(Fint, (apex, i, i2))
        end

        V = hcat(X, Y, Z)
        F = reshape(vcat(Fext, Fint), 3, :)'
        return V, F
    end
end

function calculateNormals(V::AbstractMatrix{<:Real}, F::AbstractMatrix{<:Integer})
    M = size(F, 1)
    normals = zeros(Float64, 3, M)
    areas = zeros(Float64, M)
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

        Atri = 0.5 * norm(n)
        fc = (p1 .+ p2 .+ p3) ./ 3

        if dot(n, fc .- center) < 0
            n .= -n
        end

        normals[:, i] .= n ./ norm(n)
        areas[i] = Atri
    end

    return normals, areas, sum(areas)
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

function random_fragment(obj_size::Real)
    num_points = rand(10:30)
    pts = generate_points_in_sphere(num_points)
    scale_factor = ((3 * obj_size) / (4π))^(1 / 3)
    pts .= pts .* scale_factor
    hull = qhull(pts', "Qt")
    V = pts
    F = permutedims(hull.facets, (2, 1))
    return Matrix{Float64}(V), Matrix{Int}(F)
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
            xdim = 0.5 + 0.1 * randn()
            ydim = 0.5 + 0.1 * randn()
            zdim = 0.5 + 0.1 * randn()
            siz2 = 2.1
        elseif sub_cat == 2
            xdim = 0.7 + 0.1 * randn()
            ydim = 0.7 + 0.1 * randn()
            zdim = 0.7 + 0.1 * randn()
            siz2 = 2.2
        elseif sub_cat == 3
            xdim = 1.0 + 0.1 * randn()
            ydim = 1.0 + 0.1 * randn()
            zdim = 1.0 + 0.1 * randn()
            siz2 = 2.3
        else
            xdim = 1.1; ydim = 1.5; zdim = 1.2; siz2 = 2.4
        end
        siz = 2
    else
        if sub_cat == 1
            xdim = 1.25 + 0.15 * randn()
            ydim = 1.25 + 0.15 * randn()
            zdim = 1.25 + 0.15 * randn()
            siz2 = 3.1
        elseif sub_cat == 2
            xdim = 1.5 + 0.15 * randn()
            ydim = 1.5 + 0.15 * randn()
            zdim = 1.5 + 0.15 * randn()
            siz2 = 3.2
        elseif sub_cat == 3
            xdim = 1.8 + 0.15 * randn()
            ydim = 1.8 + 0.15 * randn()
            zdim = 1.8 + 0.15 * randn()
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
            xdim = 0.1 + 0.05 * randn()
        elseif sub_cat == 2
            xdim = 0.2 + 0.05 * randn()
        elseif sub_cat == 3
            xdim = 0.3 + 0.05 * randn()
        else
            xdim = 0.1 * rand(1:5)
        end
        siz = 1
        siz2 = 1.0 + 0.1 * sub_cat
    elseif sat_size == 2
        if sub_cat == 1
            xdim = 0.6 + 0.05 * randn()
        elseif sub_cat == 2
            xdim = 0.7 + 0.05 * randn()
        elseif sub_cat == 3
            xdim = 0.8 + 0.05 * randn()
        else
            xdim = 0.9 + 0.05 * randn()
        end
        siz = 2
        siz2 = 2.0 + 0.1 * sub_cat
    else
        if sub_cat == 1
            xdim = 1.0 + 0.05 * randn()
        elseif sub_cat == 2
            xdim = 1.5 + 0.05 * randn()
        elseif sub_cat == 3
            xdim = 2.0 + 0.05 * randn()
        else
            xdim = 2.5 + 0.05 * randn()
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
            xdim = 0.1 + 0.05 * randn()
        elseif sub_cat == 2
            xdim = 0.2 + 0.05 * randn()
        elseif sub_cat == 3
            xdim = 0.3 + 0.05 * randn()
        else
            xdim = 0.1 * rand(1:5)
        end
        siz = 1
        siz2 = 1.0 + 0.1 * sub_cat
    elseif sat_size == 2
        if sub_cat == 1
            xdim = 0.6 + 0.05 * randn()
        elseif sub_cat == 2
            xdim = 0.7 + 0.05 * randn()
        elseif sub_cat == 3
            xdim = 0.8 + 0.05 * randn()
        else
            xdim = 0.9 + 0.05 * randn()
        end
        siz = 2
        siz2 = 2.0 + 0.1 * sub_cat
    else
        if sub_cat == 1
            xdim = 1.0 + 0.05 * randn()
        elseif sub_cat == 2
            xdim = 1.5 + 0.05 * randn()
        elseif sub_cat == 3
            xdim = 2.0 + 0.05 * randn()
        else
            xdim = 2.5 + 0.05 * randn()
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