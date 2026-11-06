# =============================================================================
# test_visibility.jl
#
# Pruebas del módulo Visibility con casos geométricos calculados a mano: eclipse
# (umbra, penumbra e iluminado frente a los radios analíticos de los conos), línea
# de visión (Tierra en medio, mismo lado, tangencia, extremos bajo la superficie) y
# cono de deslumbramiento solar.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_visibility.jl
#
# Es ligero: solo carga Constants y Visibility. Sol en +x ECI en todas las pruebas
# de eclipse; distancias en km.
#
# Funciones: ninguna (solo pruebas)
# =============================================================================

using Test
using LinearAlgebra
using StaticArrays

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Constants)  || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :Visibility) || include(joinpath(ROOT, "src", "Visibility.jl"))

using .Constants, .Visibility

const SUN_X = SVector(1.0, 0.0, 0.0)                 # dirección del Sol: +x ECI

@testset "Visibility" begin

    @testset "Eclipse (conos de umbra y penumbra)" begin
        # Radios analíticos de los conos a x = 7000 km detrás de la Tierra
        L_u = D_sun * R / (R_sun - R)
        L_p = D_sun * R / (R_sun + R)
        r_u = R * (1 - 7000 / L_u)                   # ≈ 6345.9 km
        r_p = R * (1 + 7000 / L_p)                   # ≈ 6411.0 km

        # En el eje de la sombra: umbra; en el lado del Sol: iluminado
        @test eclipseEarthECI(SVector(-7000.0, 0, 0), SUN_X) == (true, "umbra")
        @test eclipseEarthECI(SVector(7000.0, 0, 0), SUN_X) == (false, "lit")

        # A los lados del eje: justo dentro y fuera de cada cono
        @test eclipseEarthECI(SVector(-7000.0, r_u - 1, 0), SUN_X) == (true, "umbra")
        @test eclipseEarthECI(SVector(-7000.0, r_u + 1, 0), SUN_X) == (true, "penumbra")
        @test eclipseEarthECI(SVector(-7000.0, r_p - 1, 0), SUN_X) == (true, "penumbra")
        @test eclipseEarthECI(SVector(-7000.0, r_p + 1, 0), SUN_X) == (false, "lit")

        # Por encima del polo, perpendicular al Sol (x = 0): iluminado
        @test eclipseEarthECI(SVector(0.0, 0, 7000), SUN_X) == (false, "lit")
    end

    @testset "Línea de visión" begin
        # Tierra en medio: bloqueada; mismo lado: visible
        @test !losECI(SVector(7000.0, 0, 0), SVector(-7000.0, 0, 0))
        @test losECI(SVector(7000.0, 0, 0), SVector(7000.0, 1000, 0))

        # Segmento que pasa a 10 km por encima o por debajo de la superficie
        h = R + 10.0
        @test losECI(SVector(-8000.0, h, 0), SVector(8000.0, h, 0))
        h = R - 10.0
        @test !losECI(SVector(-8000.0, h, 0), SVector(8000.0, h, 0))

        # Los dos puntos al otro lado del mismo hemisferio (la recta corta la esfera fuera del segmento)
        @test losECI(SVector(7000.0, 0, 0), SVector(9000.0, 0, 0))

        # Extremo 50 m bajo la superficie (se proyecta) y 200 km bajo ella (sin LOS)
        @test losECI(SVector(R - 0.05, 0, 0), SVector(R + 500, 0, 0))
        @test !losECI(SVector(R - 200, 0, 0), SVector(R + 500, 0, 0))

        # Puntos coincidentes: sin LOS
        @test !losECI(SVector(7000.0, 0, 0), SVector(7000.0, 0, 0))
    end

    @testset "Cono de deslumbramiento" begin
        r_obs = SVector(7000.0, 0, 0)
        r_obj = SVector(7000.0, 1000, 0)             # visto en la dirección +y

        # Sol a 10° de la línea de visión, con cono de 15°: deslumbrado
        u_sun = SVector(sind(10), cosd(10), 0.0)
        blinded, ang = sunGlareConeECI(u_sun, r_obs, r_obj, 15.0)
        @test blinded
        @test ang ≈ 10.0 atol = 1e-9

        # Sol a 90°: no deslumbra
        blinded, ang = sunGlareConeECI(SUN_X, r_obs, r_obj, 15.0)
        @test !blinded
        @test ang ≈ 90.0 atol = 1e-9
    end
end
