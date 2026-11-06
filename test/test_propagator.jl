# =============================================================================
# test_propagator.jl
#
# Pruebas del módulo Propagator: la trayectoria propagada debe tener exactamente
# los instantes de la malla 0:Δt:T que usa scripts/main.jl para la dirección del
# Sol y las estaciones. Si el propagador devuelve un instante de más (p. ej. el
# final T cuando no es múltiplo de Δt), la fotometría lee fuera de sun_dir.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_propagator.jl
#
# Carga Types, Constants, Transformations y Propagator (DifferentialEquations
# tarda en compilar la primera vez). La órbita de prueba es una LEO circular.
#
# Funciones: propagate_leo
# =============================================================================

using Test
using Dates
using StaticArrays

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Types)           || include(joinpath(ROOT, "src", "Types.jl"))
isdefined(Main, :Constants)       || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :Transformations) || include(joinpath(ROOT, "src", "Transformations.jl"))
isdefined(Main, :Propagator)      || include(joinpath(ROOT, "src", "Propagator.jl"))

using .Types, .Constants, .Transformations, .Propagator

"""
    propagate_leo(duration_s, step_s) -> TrajectoryData

Propaga una LEO circular de prueba (a = 7000 km, i = 50°) durante `duration_s`
con salida cada `step_s`.

Por qué existe: las pruebas solo cambian la duración y el paso; así cada caso
ocupa una línea.

# Argumentos
- `duration_s::Real`: duración de la propagación (s).
- `step_s::Real`: paso de salida (s).

# Devuelve
- `TrajectoryData` de la órbita de prueba.

# Excepciones
- Las de `propagate`.
"""
function propagate_leo(duration_s::Real, step_s::Real)
    # Elementos con los mismos campos que el JSON del simulador (SMA en km, ángulos en °)
    elements = (SMA = 7000.0, ECC = 0.0, INC = 50.0, RAAN = 0.0, AOP = 0.0, TA = 0.0)
    objects = Dict(:leo => kepler_to_cartesian(elements))
    epoch = DateTime("2000-01-01T11:59:28")                    # misma época que los ejemplos
    return propagate(objects, epoch, Float64(duration_s), Float64(step_s))[:leo]
end

@testset "Propagator" begin
    @testset "Instantes de salida = malla 0:Δt:T de main.jl" begin
        # T múltiplo de Δt: 720 s con paso de 120 s -> 7 instantes (0, 120, ..., 720)
        traj = propagate_leo(720, 120)
        @test traj.t ≈ collect(0.0:120.0:720.0)

        # T no múltiplo de Δt: 900 s con paso de 120 s -> 8 instantes (0, ..., 840)
        traj = propagate_leo(900, 120)
        @test length(traj.t) == length(0.0:120.0:900.0)
        @test traj.t ≈ collect(0.0:120.0:900.0)

        # Posiciones y velocidades: una por instante
        @test length(traj.r) == length(traj.v) == length(traj.t)
    end

    @testset "Coherencia con Constants" begin
        # J2_EARTH (inclinación SSO del generador) es -C20 de la tabla del propagador
        C, _ = Propagator.load_coefs_n4()
        @test J2_EARTH == -C[3, 1]
    end
end
