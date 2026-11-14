module Integrators

# =============================================================================
# Módulo Integrators
#
# Integradores numéricos de paso fijo para ecuaciones diferenciales ordinarias
# ẋ = f(x, t). Son genéricos: el estado puede ser un número, un SVector o
# cualquier tipo con suma y producto por escalar.
#
# Por qué existe: la cinemática de la actitud (Attitude/Quaternions.jl) necesita
# un integrador de paso fijo que controlamos nosotros (subpasos y renormalización
# del cuaternión). Va en un módulo aparte para que lo reutilicen otros modelos
# (p. ej. un futuro propagador en elementos equinocciales). La órbita se sigue
# integrando con DifferentialEquations (Propagator.jl, Vern9 adaptativo).
#
# Funciones:
#   runge_kutta_4 : un paso del método de Runge-Kutta clásico de cuarto orden
#
# Dependencias: ninguna (se puede incluir antes que cualquier otro módulo).
# =============================================================================

export runge_kutta_4

"""
    runge_kutta_4(f, x, t, h) -> typeof(x)

Avanza un paso `h` la solución de `ẋ = f(x, t)` con el método de Runge-Kutta
clásico de cuarto orden:

    k1 = f(x, t)
    k2 = f(x + h/2 · k1, t + h/2)
    k3 = f(x + h/2 · k2, t + h/2)
    k4 = f(x + h · k3, t + h)
    x(t + h) ≈ x + h/6 · (k1 + 2 k2 + 2 k3 + k4)

Error local O(h⁵) y global O(h⁴): al dividir `h` entre 2 el error global baja
unas 16 veces. Si `f` solo depende de `t`, el paso es la regla de Simpson y es
exacto cuando la solución es un polinomio de grado ≤ 4.

Por qué existe: integrador de paso fijo común para los modelos que no usan
DifferentialEquations (cinemática del cuaternión en Attitude).

# Argumentos
- `f`: función `f(x, t)` que devuelve ẋ, del mismo tipo que `x`.
- `x`: estado en el instante `t`.
- `t::Real`: instante actual (s).
- `h::Real`: paso de integración (s).

# Devuelve
- El estado en `t + h`, del mismo tipo que `x`.

# Excepciones
- Ninguna.

# Ejemplo
```julia
runge_kutta_4((x, t) -> -x, 1.0, 0.0, 0.1)      # ≈ exp(-0.1)
```
"""
function runge_kutta_4(f, x, t::Real, h::Real)
    k1 = f(x, t)
    k2 = f(x + (h / 2) * k1, t + h / 2)
    k3 = f(x + (h / 2) * k2, t + h / 2)
    k4 = f(x + h * k3, t + h)
    return x + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4)
end

end # module Integrators
