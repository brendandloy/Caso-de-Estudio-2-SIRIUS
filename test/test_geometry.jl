# =============================================================================
# test_geometry.jl
#
# Pruebas del módulo Geometry: mallas (áreas frente a sus valores analíticos,
# normales sin NaN, envolvente convexa), longitud característica y ejes cuerpo,
# lectura y validación del bloque "shape", sorteo con el modelo de población
# (reproducibilidad, L_c mínima, proporciones de clases) y generate_debris.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_geometry.jl
#
# Es ligero: solo carga Types, Constants, Parameters, Transformations, Attitude
# (lectores de bloques) y Geometry. Las mallas de prueba usan dimensiones
# arbitrarias: las fórmulas no dependen de la unidad.
#
# Funciones: extents
# =============================================================================

using Test
using LinearAlgebra
using Random
using Statistics

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Types)           || include(joinpath(ROOT, "src", "Types.jl"))
isdefined(Main, :Constants)       || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :Transformations) || include(joinpath(ROOT, "src", "Transformations.jl"))
isdefined(Main, :Parameters)      || include(joinpath(ROOT, "src", "Parameters.jl"))
isdefined(Main, :JSONUtils)       || include(joinpath(ROOT, "src", "JSON_utils.jl"))
isdefined(Main, :Attitude)        || include(joinpath(ROOT, "src", "Attitude", "Attitude.jl"))
isdefined(Main, :Geometry)        || include(joinpath(ROOT, "src", "Geometry", "Geometry.jl"))

using .Parameters, .Geometry
const G = Geometry                     # acceso corto a las funciones internas

# Extensión (máx - mín) de los vértices según x, y, z
extents(V) = vec(maximum(V, dims = 1) .- minimum(V, dims = 1))

@testset "Geometry" begin

    @testset "Mallas: áreas analíticas" begin
        r = 2.0
        V, F = G.sphere_mesh(r, 20, 20)
        n, _, Atot = G.calculateNormals(V, F)
        @test size(F, 1) == 720                         # 17 bandas × 40 + 2 casquetes × 20
        @test !any(isnan, n)                            # sin triángulos degenerados en los polos
        @test isapprox(Atot, 4π * r^2; rtol = 0.03)     # poliedro inscrito: algo menos que 4πr²

        V, F = G.cube_vertices_faces(1.0, 2.0, 3.0)
        @test G.calculateNormals(V, F)[3] ≈ 2 * (2 + 6 + 3)

        R, h = 1.0, 2.0
        _, _, Atot = G.calculateNormals(G.hollow_cone(R, h)...)
        @test isapprox(Atot, π * R * hypot(R, h); rtol = 0.01)            # cono sin base
        rt = 0.4
        V, F = G.hollow_cone(R, h; top_radius = rt)
        _, _, Atot = G.calculateNormals(V, F)
        @test size(F, 1) == 100
        @test isapprox(Atot, π * (R + rt) * hypot(R - rt, h); rtol = 0.01) # tronco sin bases
        @test_throws ErrorException G.hollow_cone(R, h; top_radius = R)
    end

    @testset "Envolvente convexa y fragmento" begin
        # Tetraedro con un punto interior: 4 caras, el punto interior no aparece
        P = [0.0 0 0; 1 0 0; 0 1 0; 0 0 1; 0.1 0.1 0.1]
        F = G.convex_hull(P)
        @test size(F, 1) == 4 && !(5 in F)

        for Lc in (1e-4, 3e-4, 5e-4)                    # 0.1, 0.3 y 0.5 m en km
            V, F = G.random_fragment(Lc; rng = MersenneTwister(1))
            n, _, Atot = G.calculateNormals(V, F)
            e = extents(V)
            @test mean(e) ≈ Lc rtol = 1e-10              # L_c = (X + Y + Z)/3
            @test issorted(e; rev = true)                # ejes principales: x mayor, z menor
            @test !any(isnan, n) && Atot > 0
        end

        # Misma semilla, mismo fragmento
        @test G.random_fragment(1.0; rng = MersenneTwister(5))[1] == G.random_fragment(1.0; rng = MersenneTwister(5))[1]
    end

    @testset "Longitud característica y ejes cuerpo" begin
        # L_c de las medidas frente a (X + Y + Z)/3 de la malla generada
        for (type, dims) in (("sphere", Dict("radius_m" => 0.5)),
                             ("prism", Dict("x_m" => 3.0, "y_m" => 2.0, "z_m" => 1.0)),
                             ("plate", Dict("x_m" => 2.0, "y_m" => 1.0, "thickness_m" => 0.02)),
                             ("cylinder", Dict("radius_m" => 1.0, "height_m" => 6.0)),
                             ("cone", Dict("radius_m" => 0.5, "height_m" => 1.5)),
                             ("trunk", Dict("radius_m" => 1.0, "top_radius_m" => 0.5, "height_m" => 1.0)))
            e = G.SHAPES[type]
            V, _ = e.mesh(dims, MersenneTwister(1))
            @test isapprox(e.lc(dims), mean(extents(V)); rtol = 0.02)   # malla poligonal: casi exacto
        end

        # Ejes: cilindro con eje z centrado, cono con la base en z = 0 y el vértice en +z
        V, _ = G.cylinder_mesh(1.0, 6.0)
        @test extrema(V[:, 3]) == (-3.0, 3.0)
        V, _ = G.hollow_cone(0.5, 1.5)
        @test extrema(V[:, 3]) == (0.0, 1.5)
    end

    @testset "Bloque shape: lectura y validación" begin
        ok = Dict("type" => "cylinder", "radius_m" => 1.5, "height_m" => 8.0, "seed" => 1)
        s = parse_shape(ok, "deb")
        @test s.type == "cylinder" && s.seed == 1 && char_length(s) ≈ (6 + 8) / 3

        @test_throws ErrorException shape_block(Dict("SMA" => 7000), "deb")                       # sin shape
        @test_throws ErrorException parse_shape(Dict("type" => "cube", "x_m" => 1), "deb")         # forma desconocida
        @test_throws ErrorException parse_shape(Dict("type" => "cylinder", "radius_m" => 1.0), "deb")   # medidas incompletas
        @test_throws ErrorException parse_shape(merge(ok, Dict("height_m" => -1.0)), "deb")        # medida negativa
        @test_throws ErrorException parse_shape(merge(ok, Dict("radio_m" => 1.0)), "deb")          # clave desconocida
        @test_throws ErrorException parse_shape(Dict("type" => "trunk", "radius_m" => 1.0,
            "top_radius_m" => 1.0, "height_m" => 1.0, "seed" => 1), "deb")                        # r superior = R
        @test_throws ErrorException parse_shape(Dict("type" => "plate", "x_m" => 1.0,
            "y_m" => 0.5, "thickness_m" => 0.8, "seed" => 1), "deb")                              # espesor no mínimo
        @test_throws ErrorException resolve_shape(Dict("type" => "cylinder", "radius_m" => 1.0), "deb"; seed = 1)

        # Medidas concretas por debajo del mínimo: se aceptan con un aviso
        small = Dict("type" => "sphere", "radius_m" => 0.02, "seed" => 1)
        @test_logs (:warn, r"below SHAPE_LC_MIN_M") parse_shape(small, "deb")
    end

    @testset "Sorteo con el modelo de población" begin
        # Reproducible: la misma semilla da el mismo bloque
        @test resolve_shape(Dict("type" => "random"), "deb"; seed = 9) ==
              resolve_shape(Dict("type" => "random"), "deb"; seed = 9)

        # Forma elegida sin medidas: se sortean todas y respetan las proporciones maestras
        for seed in 1:200
            b = resolve_shape(Dict("type" => "cylinder"), "deb"; seed = seed)
            hr = b["height_m"] / b["radius_m"]
            @test CYLINDER_HR_RANGE[1] - 1e-4 <= hr <= CYLINDER_HR_RANGE[2] + 1e-4
        end

        # Prisma sorteado: z es la arista mayor (dirección longitudinal) y x ≥ y
        for seed in 1:200
            b = resolve_shape(Dict("type" => "prism"), "deb"; seed = seed)
            @test b["z_m"] >= b["x_m"] >= b["y_m"]
        end

        # 4000 objetivos "random": todos con L_c ≥ L_min y clases en proporción a la tabla 3.1
        N = 4000
        blocks = [resolve_shape(Dict("type" => "random"), "deb"; seed = s) for s in 1:N]
        lcs = [char_length(parse_shape(b, "deb")) for b in blocks]
        @test all(>=(SHAPE_LC_MIN_M * (1 - 1e-5)), lcs)
        c = SHAPE_CLASS_COUNTS
        p_frag = c.fragment / (c.fragment + c.payload + c.rocket_body + c.mission_related)
        @test isapprox(mean(b["class"] == "fragment" for b in blocks), p_frag; atol = 0.03)

        # Fragmentos: N(>L) ∝ L^-1.6 entre 10 cm y 1 m (fracción por debajo de 0.3 m)
        frag = [b["lc_m"] for b in blocks if b["type"] == "fragment"]
        β = FRAGMENT_LC_BETAS[3]
        p03 = (0.1^-β - 0.3^-β) / (0.1^-β - 1.0^-β)
        @test isapprox(mean(frag .< 0.3), p03; atol = 0.04)

        # Con L_min = 1 cm, la gran mayoría serían fragmentos (≈ 92 %)
        w = G.class_weights(0.01)
        @test w.fragment / sum(values(w)) > 0.9
    end

    @testset "generate_debris: todas las formas" begin
        for type in vcat(sort!(collect(keys(G.SHAPES))), "random"), seed in 1:5
            s = parse_shape(resolve_shape(Dict("type" => type), "deb"; seed = seed), "deb")
            d = generate_debris(s)
            @test !any(isnan, d[:normals]) && d[:total_area] > 0
            @test size(d[:normals], 2) == length(d[:areas])
            @test d[:C_spec] + d[:C_diff] <= 1 + 1e-12
            @test generate_debris(s)[:areas] == d[:areas]              # misma semilla, mismo objeto
        end
    end
end
