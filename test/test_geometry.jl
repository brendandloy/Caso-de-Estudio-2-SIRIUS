# =============================================================================
# test_geometry.jl
#
# Pruebas de las mallas de Geometry.jl: áreas frente a sus valores analíticos,
# longitud característica del fragmento, ausencia de normales NaN y de
# dimensiones negativas, y generación de todas las formas.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_geometry.jl
#
# Solo carga Geometry.jl (necesita QHull, que usa scipy a través de PyCall).
# Las dimensiones de las mallas son arbitrarias: las fórmulas no dependen de la unidad.
#
# Funciones: char_length
# =============================================================================

using Test
using LinearAlgebra
using Random
using Statistics

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))
isdefined(Main, :Geometry) || include(joinpath(ROOT, "src", "Geometry.jl"))
using .Geometry
const G = Geometry                     # acceso corto a las funciones internas

"""
    char_length(V) -> Float64

Longitud característica `(X + Y + Z)/3` de una nube de puntos, con X, Y, Z sus
dimensiones en los ejes principales (la misma definición que usa `random_fragment`).

# Argumentos
- `V::AbstractMatrix`: puntos (N × 3).

# Devuelve
- `Float64`.

# Excepciones
- Ninguna.
"""
function char_length(V::AbstractMatrix)
    P = V .- mean(V, dims = 1)
    proj = P * eigen(Symmetric(P' * P)).vectors
    return mean(vec(maximum(proj, dims = 1) .- minimum(proj, dims = 1)))
end

Random.seed!(1)                         # resultados reproducibles

@testset "Geometry" begin

    @testset "Esfera" begin
        r = 2.0
        V, F = G.sphere_mesh(r, 20, 20)
        n, A, Atot = G.calculateNormals(V, F)
        @test size(F, 1) == 720                         # 17 bandas × 40 + 2 casquetes × 20
        @test !any(isnan, n)                            # sin triángulos degenerados en los polos
        @test isapprox(Atot, 4π * r^2; rtol = 0.03)     # poliedro inscrito: algo menos que 4πr²
        @test all(dot(n[:, k], vec(mean(V[F[k, :], :], dims = 1))) > 0 for k in axes(F, 1))   # normales hacia fuera
    end

    @testset "Cono y tronco (sin caras duplicadas)" begin
        R, h = 1.0, 2.0
        V, F = G.hollow_cone(R, h, false)
        _, _, Atot = G.calculateNormals(V, F)
        @test size(F, 1) == 50
        @test isapprox(Atot, π * R * sqrt(R^2 + h^2); rtol = 0.01)      # superficie lateral

        V, F = G.hollow_cone(R, h, true)
        _, _, Atot = G.calculateNormals(V, F)
        r = norm(V[51, 1:2])                                            # radio superior sorteado
        @test size(F, 1) == 100
        @test isapprox(Atot, π * (R + r) * sqrt((R - r)^2 + h^2); rtol = 0.01)
    end

    @testset "Fragmento: longitud característica" begin
        for Lc in (1e-4, 3e-4, 5e-4)                     # 0.1, 0.3 y 0.5 m en km
            V, F = G.random_fragment(Lc)
            n, _, Atot = G.calculateNormals(V, F)
            @test char_length(V) ≈ Lc rtol = 1e-10
            @test maximum(norm.(eachrow(V))) < 2Lc       # nada de radios de decenas de metros
            @test !any(isnan, n) && Atot > 0
        end
    end

    @testset "Triángulos degenerados" begin
        V = [0.0 0 0; 1 0 0; 0 1 0; 2 0 0]
        F = [1 2 3; 1 2 4]                               # el segundo tiene los tres vértices alineados
        n, A, _ = @test_logs (:warn,) G.calculateNormals(V, F)
        @test size(n, 2) == 1 && length(A) == 1
    end

    @testset "Dimensiones positivas" begin
        # 2000 sorteos por categoría; antes ≈ 2–3 % de los conos pequeños salían con dimensión negativa
        draws(f) = [x for sat in 1:3, sub in 1:4, _ in 1:2000 for x in f(sat, sub)[1:3]]
        @test all(>(0), draws(G._dims_prism))
        @test all(>(0), draws(G._dims_cone))
        @test all(>(0), draws(G._dims_trunk))
    end

    @testset "generate_debris: todas las formas" begin
        for shape in 1:7, sat in 1:3, sub in 1:4
            d = generate_debris(sat, sub; shape_type = shape)
            @test !any(isnan, d[:normals]) && d[:total_area] > 0
            @test size(d[:normals], 2) == length(d[:areas])
        end
    end
end
