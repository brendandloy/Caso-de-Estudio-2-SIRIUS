# =============================================================================
# test_transformations.jl
#
# Pruebas del módulo Transformations con casos calculados a mano: elementos
# keplerianos -> estado cartesiano (perigeo, apogeo, energía, momento angular,
# inclinación y nodo), ángulo de rotación terrestre (ERA), rotación ECEF -> ECI
# y coordenadas geodésicas -> ECEF (WGS84).
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_transformations.jl
#
# Es ligero: solo carga Types, Constants y Transformations. Unidades: km, km/s, grados.
#
# Funciones: elements (solo construye la entrada de kepler_to_cartesian)
# =============================================================================

using Test
using LinearAlgebra
using StaticArrays

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Types)           || include(joinpath(ROOT, "src", "Types.jl"))
isdefined(Main, :Constants)       || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :Transformations) || include(joinpath(ROOT, "src", "Transformations.jl"))

using .Types, .Constants, .Transformations

# Elementos con los mismos campos que el JSON del simulador (SMA en km, ángulos en °)
elements(a, e, i, Ω, ω, ν) = (SMA = a, ECC = e, INC = i, RAAN = Ω, AOP = ω, TA = ν)

@testset "Transformations" begin

    @testset "Kepler -> cartesianas" begin
        # Circular ecuatorial: r sobre +x, v = √(μ/a) sobre +y; con i = 90°, v sobre +z
        s = kepler_to_cartesian(elements(7000.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        @test s.r ≈ SVector(7000.0, 0, 0)
        @test s.v ≈ SVector(0, sqrt(μ / 7000), 0)
        @test kepler_to_cartesian(elements(7000.0, 0.0, 90.0, 0.0, 0.0, 0.0)).v ≈ SVector(0, 0, sqrt(μ / 7000))

        # Elíptica: perigeo a(1 - e) con v_p = √(μ/a · (1+e)/(1-e)); apogeo a(1 + e) en -x
        a, e = 10000.0, 0.3
        s = kepler_to_cartesian(elements(a, e, 0.0, 0.0, 0.0, 0.0))
        @test s.r ≈ SVector(a * (1 - e), 0, 0)
        @test norm(s.v) ≈ sqrt(μ / a * (1 + e) / (1 - e))
        @test kepler_to_cartesian(elements(a, e, 0.0, 0.0, 0.0, 180.0)).r ≈ SVector(-a * (1 + e), 0, 0)

        # Órbita general: energía -μ/2a, |h| = √(μ a (1 - e²)) y ĥ = (sin i sin Ω, -sin i cos Ω, cos i)
        a, e, i, Ω = 12000.0, 0.25, 63.4, 40.0
        s = kepler_to_cartesian(elements(a, e, i, Ω, 270.0, 123.0))
        h = cross(s.r, s.v)
        @test norm(s.v)^2 / 2 - μ / norm(s.r) ≈ -μ / (2a)
        @test norm(h) ≈ sqrt(μ * a * (1 - e^2))
        @test acosd(h[3] / norm(h)) ≈ i
        @test atand(h[1], -h[2]) ≈ Ω
    end

    @testset "Rotación terrestre y sistemas de referencia" begin
        # ERA en J2000: 2π · 0.7790572732640 rad = 280.46061837504°; avanza 2π · 1.00273781191135448 por día UT
        @test rad2deg(theta0_from_utc(jd_J2000)) ≈ 280.46061837504
        @test mod(theta0_from_utc(jd_J2000 + 1) - theta0_from_utc(jd_J2000), 2π) ≈ 2π * 0.00273781191135448

        # ECEF -> ECI: con la Tierra girada 90°, el eje x ECEF apunta a +y ECI
        @test rotate_z(SVector(1.0, 0, 0), π / 2) ≈ SVector(0, 1.0, 0)

        # WGS84: ecuador a 0° y 1 km de altura -> a + 1; polo -> b = a (1 - f)
        @test llh_to_ecef(0.0, 0.0, 1.0) ≈ SVector(WGS84_A + 1, 0, 0)
        @test llh_to_ecef(π / 2, 0.0, 0.0) ≈ SVector(0, 0, WGS84_A * (1 - WGS84_f)) atol = 1e-9
    end
end
