# =============================================================================
# runtests.jl
#
# Ejecuta todas las pruebas de SIRIUS. Cada archivo se carga en su propio módulo
# para que sus constantes (ROOT, trayectorias de prueba...) y los módulos de src/
# que incluye no choquen con los de los demás.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. test/runtests.jl
#   o, en el REPL con el entorno activo:  ] test
#
# Cada archivo se puede ejecutar también por separado (ver su cabecera).
#
# Funciones: ninguna (solo pruebas)
# =============================================================================

using Test

# Archivos de prueba, de los módulos más básicos a los que dependen de más módulos
const TEST_FILES = [
    "test_transformations.jl",
    "test_propagator.jl",
    "test_visibility.jl",
    "test_attitude.jl",
    "test_geometry.jl",
    "test_photometry.jl",
    "test_camera.jl",
    "test_scenario.jl",
]

@testset "SIRIUS" begin
    for file in TEST_FILES
        # Módulo propio por archivo, p. ej. Test_test_attitude
        name = Symbol("Test_", splitext(file)[1])
        path = joinpath(@__DIR__, file)
        @eval module $name
            include($path)
        end
    end
end
