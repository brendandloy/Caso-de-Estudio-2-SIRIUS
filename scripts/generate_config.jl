# =============================================================================
# generate_config.jl
#
# Script de entrada del generador de escenarios: lee una spec corta y escribe
# el JSON completo que consume el simulador.
#
# Uso:
#   julia --project=. scripts/generate_config.jl [<spec.json> [<salida.json>]]
#
# Por defecto: config/specs/ejemplo_leo.json -> config/generated/ejemplo_leo.json
#
# Es ligero a propósito: solo carga Types, Constants, Transformations, Attitude y Scenario
# (no GLMakie ni DifferentialEquations). Scenario hace `using ..Constants`,
# `using ..Transformations` y `using ..Attitude`, que aquí resuelven a los módulos de Main.
# Transformations aporta theta0_from_utc, la misma rotación terrestre que usa el simulador,
# y Attitude valida los bloques "attitude" con el mismo código que el simulador.
#
# Funciones: run_generator
# =============================================================================

# Raíz del repositorio (este script vive en scripts/)
const ROOT = normpath(joinpath(@__DIR__, ".."))

# Orden obligatorio: Transformations depende de Types y Constants
isdefined(Main, :Types)           || include(joinpath(ROOT, "src", "Types.jl"))
isdefined(Main, :Constants)       || include(joinpath(ROOT, "src", "Constants.jl"))
isdefined(Main, :Transformations) || include(joinpath(ROOT, "src", "Transformations.jl"))
isdefined(Main, :Attitude)        || include(joinpath(ROOT, "src", "Attitude", "Attitude.jl"))
include(joinpath(ROOT, "src", "Scenario", "Scenario.jl"))

"""
    run_generator(args=ARGS) -> String

Genera el JSON de escenario a partir de los argumentos de línea de comandos.

Por qué existe: separa la lógica de los argumentos del `include`, de modo que se
puede llamar también desde el REPL sin lanzar el script entero.

# Argumentos
- `args`: `[spec.json [salida.json]]`. Sin argumentos usa la spec de ejemplo;
  sin salida escribe en `config/generated/<nombre de la spec>.json`.

# Devuelve
- `String`: la ruta del JSON generado.

# Excepciones
- `ErrorException`: si la spec no existe o no es válida (ver `Scenario.generate_config_file`).

# Ejemplo
```julia
run_generator(["config/specs/ejemplo_sso_geo.json"])
# escribe config/generated/ejemplo_sso_geo.json

run_generator(["config/specs/ejemplo_sso_geo.json", "config/mi_escenario.json"])
# escribe config/mi_escenario.json
```
"""
function run_generator(args = ARGS)
    # Spec de entrada: la indicada o la de ejemplo
    spec_path = isempty(args) ?
        joinpath(ROOT, "config", "specs", "ejemplo_leo.json") :
        abspath(args[1])

    # Salida: la indicada o config/generated/ con el nombre de la spec
    out_path = length(args) >= 2 ?
        abspath(args[2]) :
        joinpath(ROOT, "config", "generated", splitext(basename(spec_path))[1] * ".json")

    return Scenario.generate_config_file(spec_path, out_path)
end

# Solo se ejecuta al lanzar el script directamente, no al hacer include desde el REPL
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_generator()
end