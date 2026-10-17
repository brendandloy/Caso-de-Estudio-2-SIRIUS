# Uso:
#   julia --project=. scripts/generate_config.jl [<spec.json> [<salida.json>]]
#
# Por defecto: config/specs/ejemplo_leo.json -> config/generated/ejemplo_leo.json
#
# Es ligero a propósito: solo carga Constants y Scenario (no GLMakie ni DifferentialEquations).
# Scenario hace `using ..Constants`, que aquí resuelve a Main.Constants.

const ROOT = normpath(joinpath(@__DIR__, ".."))

isdefined(Main, :Constants) || include(joinpath(ROOT, "src", "Constants.jl"))
include(joinpath(ROOT, "src", "Scenario", "Scenario.jl"))

function run_generator(args = ARGS)
    spec_path = isempty(args) ?
        joinpath(ROOT, "config", "specs", "ejemplo_leo.json") :
        abspath(args[1])

    out_path = length(args) >= 2 ?
        abspath(args[2]) :
        joinpath(ROOT, "config", "generated", splitext(basename(spec_path))[1] * ".json")

    return Scenario.generate_config_file(spec_path, out_path)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    run_generator()
end
