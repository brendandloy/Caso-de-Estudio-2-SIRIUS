# =============================================================================
# test_scenario.jl
#
# Pruebas del generador de escenarios (módulo Scenario): reproducibilidad (misma
# spec y semilla -> mismo JSON; otra semilla -> otro escenario), ida y vuelta
# elementos -> r, v -> elementos, misma convención que el simulador
# (kepler_to_cartesian) e inclinación de una órbita heliosíncrona.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_scenario.jl
#
# Carga Types, Constants, Parameters, JSON_utils, Transformations, Attitude,
# Geometry y Scenario. El generador imprime las formas sorteadas; es normal.
#
# Funciones: spec
# =============================================================================

using Test
using LinearAlgebra

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Types)           || include(joinpath(ROOT, "src", "Types.jl"))
isdefined(Main, :Constants)       || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :Parameters)      || include(joinpath(ROOT, "src", "Parameters.jl"))
isdefined(Main, :JSONUtils)       || include(joinpath(ROOT, "src", "JSON_utils.jl"))
isdefined(Main, :Transformations) || include(joinpath(ROOT, "src", "Transformations.jl"))
isdefined(Main, :Attitude)        || include(joinpath(ROOT, "src", "Attitude", "Attitude.jl"))
isdefined(Main, :Geometry)        || include(joinpath(ROOT, "src", "Geometry", "Geometry.jl"))
isdefined(Main, :Scenario)        || include(joinpath(ROOT, "src", "Scenario", "Scenario.jl"))

using .Constants, .Transformations, .Scenario
const S = Scenario                     # acceso corto a las funciones internas

# Spec mínima: un objetivo con órbita, actitud y forma aleatorias y una cámara en LEO
spec(seed) = Dict{String,Any}(
    "seed" => seed,
    "objectives" => Dict{String,Any}("deb" => Dict{String,Any}(
        "orbit"    => Dict{String,Any}("type" => "random_LEO"),
        "attitude" => Dict{String,Any}("type" => "random"),
        "shape"    => Dict{String,Any}("type" => "random"))),
    "observers" => Dict{String,Any}("cam" => Dict{String,Any}(
        "orbit"    => Dict{String,Any}("type" => "LEO", "altitude_km" => 600, "INC" => 97.0, "RAAN" => 0.0, "TA" => 0.0),
        "attitude" => Dict{String,Any}("type" => "fixed", "frame" => "ECI", "q" => [1, 0, 0, 0]))))

@testset "Scenario" begin

    @testset "Reproducibilidad" begin
        cfg = generate_config(spec(7))
        @test generate_config(spec(7)) == cfg                    # misma semilla: mismo escenario
        other = generate_config(spec(8))
        @test other["objectives"]["deb"]["SMA"] != cfg["objectives"]["deb"]["SMA"]   # otra semilla: otra órbita
        @test other["objectives"]["deb"]["shape"]["seed"] != cfg["objectives"]["deb"]["shape"]["seed"]
    end

    @testset "Elementos <-> r, v" begin
        el = Dict{String,Any}("SMA" => 12000.0, "ECC" => 0.25, "INC" => 63.4,
                              "RAAN" => 40.0, "AOP" => 270.0, "TA" => 123.0)
        r, v = S.elements_to_rv(el, μ)
        back = S.cartesian_to_kepler(r, v, μ)
        for k in keys(el)
            @test back[k] ≈ el[k]
        end

        # Misma convención que el simulador: el generador y kepler_to_cartesian dan el mismo estado
        s = kepler_to_cartesian((SMA = 12000.0, ECC = 0.25, INC = 63.4, RAAN = 40.0, AOP = 270.0, TA = 123.0))
        @test s.r ≈ r && s.v ≈ v
    end

    @testset "Órbita heliosíncrona" begin
        # A 700 km, la precesión nodal por J2 sigue al Sol medio con i ≈ 98,19° (circular)
        @test S.sso_inclination_deg(R + 700, 0.0, (mu = μ, R = R)) ≈ 98.19 atol = 0.01
    end
end
