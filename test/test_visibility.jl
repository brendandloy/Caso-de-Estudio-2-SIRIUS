# =============================================================================
# test_visibility.jl
#
# Pruebas del módulo Visibility con casos geométricos calculados a mano: eclipse
# (umbra, penumbra e iluminado frente a los radios analíticos de los conos), línea
# de visión (Tierra en medio, mismo lado, tangencia, extremos bajo la superficie) y
# ángulo entre el Sol y el eje óptico de la cámara (exclusión solar).
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
        L_u = D_SUN * R_EARTH / (R_SUN - R_EARTH)
        L_p = D_SUN * R_EARTH / (R_SUN + R_EARTH)
        r_u = R_EARTH * (1 - 7000 / L_u)                   # ≈ 6345.9 km
        r_p = R_EARTH * (1 + 7000 / L_p)                   # ≈ 6411.0 km

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
        h = R_EARTH + 10.0
        @test losECI(SVector(-8000.0, h, 0), SVector(8000.0, h, 0))
        h = R_EARTH - 10.0
        @test !losECI(SVector(-8000.0, h, 0), SVector(8000.0, h, 0))

        # Los dos puntos al otro lado del mismo hemisferio (la recta corta la esfera fuera del segmento)
        @test losECI(SVector(7000.0, 0, 0), SVector(9000.0, 0, 0))

        # Extremo 50 m bajo la superficie (se proyecta) y 200 km bajo ella (sin LOS)
        @test losECI(SVector(R_EARTH - 0.05, 0, 0), SVector(R_EARTH + 500, 0, 0))
        @test !losECI(SVector(R_EARTH - 200, 0, 0), SVector(R_EARTH + 500, 0, 0))

        # Puntos coincidentes: sin LOS
        @test !losECI(SVector(7000.0, 0, 0), SVector(7000.0, 0, 0))
    end

    @testset "Ángulo Sol-boresight" begin
        # C_ci = I: el eje óptico es +z ECI; Sol en +z, +x y -z
        Id = Matrix{Float64}(I, 3, 3)
        @test sun_boresight_angle_deg(Id, SVector(0.0, 0, 1)) ≈ 0.0 atol = 1e-9
        @test sun_boresight_angle_deg(Id, SUN_X) ≈ 90.0 atol = 1e-9
        @test sun_boresight_angle_deg(Id, SVector(0.0, 0, -1)) ≈ 180.0 atol = 1e-9

        # Cámara girada con el eje óptico en +y ECI (filas x_c = +x, y_c = -z, z_c = +y; det = +1)
        # y el Sol a 20° de +y en el plano xy
        C_y = [1.0 0 0; 0 0 -1; 0 1 0]
        @test det(C_y) ≈ 1.0
        @test sun_boresight_angle_deg(C_y, SVector(sind(20), cosd(20), 0.0)) ≈ 20.0 atol = 1e-9
    end
end
