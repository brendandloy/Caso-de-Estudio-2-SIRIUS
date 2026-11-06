# =============================================================================
# main.jl
#
# Script principal de SIRIUS: ejecuta la simulación completa a partir de una spec.
#
# Uso:
#   julia --project=. scripts/main.jl <spec.json> [<output_path>]
#   (desde VS Code: scripts/run_simulation_vscode.jl, cambiando la variable `spec`)
#
# Pasos (ver ?run_spec y ?run_simulation, en src/Simulation.jl):
#   1. Genera el JSON completo del escenario en config/generated/<nombre de la spec>.json
#      (si ya existe, se sobrescribe; guarda las semillas usadas, así que se puede repetir).
#   2. Simula con ese JSON: propagación -> forma y actitud de los objetivos -> fotometría ->
#      actitud de los observadores y orientación de su cámara -> cámara -> CSV.
# Sin <output_path>, los resultados se escriben en output/<nombre de la spec>/.
#
# Al hacer include (VS Code, REPL, benchmark.jl, MonteCarlo.jl) solo carga el paquete y deja
# disponibles run_spec y run_simulation; la simulación se lanza al ejecutarlo desde la terminal.
#
# Funciones: ninguna (están en src/Simulation.jl).
# =============================================================================

push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using Project
using Project.Simulation                     # run_spec, run_simulation
using Project.Environment: parse_arguments   # argumentos de la línea de comandos

# Solo se ejecuta al lanzar el script directamente, no al hacer include (VS Code, MonteCarlo)
if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    args = parse_arguments()
    run_spec(args.spec_path, args.output_path)
end
