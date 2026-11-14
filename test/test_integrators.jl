# =============================================================================
# test_integrators.jl
#
# Pruebas del módulo Integrators con problemas de solución analítica conocida:
# exactitud con polinomios (regla de Simpson), orden 4 de convergencia y un
# estado vectorial (oscilador armónico).
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/test_integrators.jl
#
# Es ligero: solo carga Integrators.
#
# Funciones: integrate_fixed
# =============================================================================

using Test
using StaticArrays

# Raíz del repositorio (este script vive en test/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Integrators) || include(joinpath(ROOT, "src", "Integrators.jl"))

using .Integrators

"""
    integrate_fixed(f, x0, T, N) -> typeof(x0)

Integra `ẋ = f(x, t)` de 0 a `T` con `N` pasos iguales de `runge_kutta_4`.

# Argumentos
- `f`: función `f(x, t)`. `x0`: estado inicial. `T::Real`: tiempo final. `N::Int`: número de pasos.

# Devuelve
- El estado en `T`.

# Excepciones
- Ninguna.
"""
function integrate_fixed(f, x0, T, N)
    x, h = x0, T / N
    for k in 0:N-1
        x = runge_kutta_4(f, x, k * h, h)
    end
    return x
end

@testset "Integrators" begin

    @testset "Exacto con polinomios" begin
        # ẋ = 4t³ -> x = t⁴: un solo paso de h = 2 es la regla de Simpson, exacta
        @test runge_kutta_4((x, t) -> 4t^3, 0.0, 0.0, 2.0) ≈ 16.0 atol = 1e-12
        # Derivada constante: x = x0 + c h
        @test runge_kutta_4((x, t) -> 3.0, 1.0, 0.0, 0.5) ≈ 2.5 atol = 1e-15
    end

    @testset "Orden 4" begin
        # ẋ = -x, x(0) = 1 -> x(1) = e⁻¹. Al dividir h entre 2, el error baja ≈ 2⁴ = 16
        f(x, t) = -x
        e10 = abs(integrate_fixed(f, 1.0, 1.0, 10) - exp(-1))
        e20 = abs(integrate_fixed(f, 1.0, 1.0, 20) - exp(-1))
        @test 15 < e10 / e20 < 18
        @test e10 < 1e-6
    end

    @testset "Estado vectorial" begin
        # Oscilador armónico ẍ = -x: tras un periodo (2π) vuelve a [1, 0]
        f(x, t) = SVector(x[2], -x[1])
        xT = integrate_fixed(f, SVector(1.0, 0.0), 2π, 1000)
        @test xT isa SVector{2,Float64}
        @test maximum(abs.(xT - SVector(1.0, 0.0))) < 1e-9
    end
end
