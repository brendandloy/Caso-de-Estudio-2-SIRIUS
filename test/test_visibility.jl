# =============================================================================
# test_visibility.jl
#
# Pruebas del módulo Visibility con casos geométricos calculados a mano: función de
# sombra (umbra, penumbra e iluminado frente a los radios analíticos de los conos, y
# el limbo de la Tierra en el centro del Sol, ν ≈ 1/2), línea
# de visión (Tierra en medio, mismo lado, atmósfera, extremos bajo la superficie),
# ángulo entre el Sol y el eje óptico de la cámara (exclusión solar) y la Tierra vista
# desde el observador (radio con la atmósfera, semiángulo del disco, choque de un rayo).
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_visibility.jl
#
# Es ligero: solo carga Constants, Parameters y Visibility. Sol en +x ECI en todas las pruebas
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
isdefined(Main, :Parameters) || include(joinpath(ROOT, "src", "Parameters.jl"))
isdefined(Main, :Visibility) || include(joinpath(ROOT, "src", "Visibility.jl"))

using .Constants, .Parameters, .Visibility

const SUN_X = SVector(1.0, 0.0, 0.0)                 # dirección del Sol: +x ECI

@testset "Visibility" begin

    @testset "Función de sombra (conos de umbra y penumbra)" begin
        # Radios analíticos de los conos a x = 7000 km detrás de la Tierra
        L_u = D_SUN * R_EARTH / (R_SUN - R_EARTH)
        L_p = D_SUN * R_EARTH / (R_SUN + R_EARTH)
        r_u = R_EARTH * (1 - 7000 / L_u)                   # ≈ 6345.9 km
        r_p = R_EARTH * (1 + 7000 / L_p)                   # ≈ 6411.0 km

        # En el eje de la sombra: umbra (ν = 0); en el lado del Sol: iluminado (ν = 1)
        @test shadow_function(SVector(-7000.0, 0, 0), SUN_X) == 0.0
        @test shadow_function(SVector(7000.0, 0, 0), SUN_X) == 1.0

        # A los lados del eje, justo dentro y fuera de cada cono: el borde de la umbra es la
        # tangencia interior de los discos (c = b - a) y el de la penumbra la exterior (c = a + b)
        @test shadow_function(SVector(-7000.0, r_u - 1, 0), SUN_X) == 0.0
        @test 0.0 < shadow_function(SVector(-7000.0, r_u + 1, 0), SUN_X) < 1.0
        @test 0.0 < shadow_function(SVector(-7000.0, r_p - 1, 0), SUN_X) < 1.0
        @test shadow_function(SVector(-7000.0, r_p + 1, 0), SUN_X) == 1.0

        # Por encima del polo, perpendicular al Sol (x = 0): iluminado
        @test shadow_function(SVector(0.0, 0, 7000), SUN_X) == 1.0

        # Limbo de la Tierra sobre el centro del Sol (c = b): se tapa casi la mitad del disco
        # solar, ν ≈ 1/2 (el limbo es casi recto a la escala del Sol, a/b ≈ 0,004)
        b = asin(R_EARTH / 7000)
        @test shadow_function(7000 * SVector(-cos(b), sin(b), 0), SUN_X) ≈ 0.5 atol = 0.01
    end

    @testset "Línea de visión" begin
        # Tierra en medio: bloqueada; mismo lado: visible
        @test !losECI(SVector(7000.0, 0, 0), SVector(-7000.0, 0, 0))
        @test losECI(SVector(7000.0, 0, 0), SVector(7000.0, 1000, 0))

        # Segmento que pasa 10 km por encima de la atmósfera (visible), 10 km por debajo de su
        # límite (tapado por la atmósfera) y 10 km bajo la superficie (tapado)
        R_atm = R_EARTH + EARTH_ATMOSPHERE_ALTITUDE_KM
        for (h, has_los) in ((R_atm + 10.0, true), (R_atm - 10.0, false), (R_EARTH - 10.0, false))
            @test losECI(SVector(-8000.0, h, 0), SVector(8000.0, h, 0)) == has_los
        end

        # Los dos puntos al otro lado del mismo hemisferio (la recta corta la esfera fuera del segmento)
        @test losECI(SVector(7000.0, 0, 0), SVector(9000.0, 0, 0))

        # Extremo 50 m bajo la superficie (se proyecta) y 200 km bajo ella (sin LOS)
        @test losECI(SVector(R_EARTH - 0.05, 0, 0), SVector(R_EARTH + 500, 0, 0))
        @test !losECI(SVector(R_EARTH - 200, 0, 0), SVector(R_EARTH + 500, 0, 0))

        # Puntos coincidentes: sin LOS
        @test !losECI(SVector(7000.0, 0, 0), SVector(7000.0, 0, 0))
    end

    @testset "Tierra vista desde el observador" begin
        # Radio que tapa: con la atmósfera si los dos puntos están por encima de ella; si no, la superficie
        R_atm = R_EARTH + EARTH_ATMOSPHERE_ALTITUDE_KM
        @test earth_occulting_radius(SVector(7000.0, 0, 0)) == R_atm
        @test earth_occulting_radius(SVector(7000.0, 0, 0), SVector(R_EARTH + 50, 0, 0)) == R_EARTH

        # A 2R del centro, el disco de la Tierra tiene un semiángulo de 30° (sin ρ = 1/2)
        @test earth_angular_radius_deg(SVector(2R_atm, 0, 0)) ≈ 30.0

        # Rayos desde x = 10000 km: hacia el centro choca en x = R; hacia fuera y en perpendicular, no
        r_obs = SVector(10000.0, 0, 0)
        @test earth_ray_hit(r_obs, SVector(-1.0, 0, 0)) ≈ SVector(R_atm, 0, 0)
        @test earth_ray_hit(r_obs, SVector(1.0, 0, 0)) === nothing
        @test !is_earth_occulted(r_obs, SVector(0.0, 1, 0))

        # Direcciones a ρ⊕ ∓ 0,01° del nadir: dentro del disco (tapada) y fuera (visible)
        ρ = earth_angular_radius_deg(r_obs)
        @test is_earth_occulted(r_obs, SVector(-cosd(ρ - 0.01), sind(ρ - 0.01), 0))
        @test !is_earth_occulted(r_obs, SVector(-cosd(ρ + 0.01), sind(ρ + 0.01), 0))
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
