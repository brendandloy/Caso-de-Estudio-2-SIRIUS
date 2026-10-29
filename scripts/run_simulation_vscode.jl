# =============================================================================
# run_simulation_vscode.jl
#
# Lanzador para VS Code: ejecuta la simulación completa (spec -> JSON completo en
# config/generated/ -> simulación) con la spec indicada abajo.
#
# Uso: cambia `spec` por la spec que quieras simular y ejecuta este archivo
# (botón ▷ de VS Code o "Julia: Execute active File in REPL").
#
# Activa el entorno del repositorio si no lo está, así funciona aunque Julia se
# lance sin --project=. (las dependencias son las de Project.toml/Manifest.toml).
# =============================================================================

# Entorno del repositorio (carpeta padre de scripts/)
import Pkg
repo_root = dirname(@__DIR__)
if Base.active_project() != joinpath(repo_root, "Project.toml")
    Pkg.activate(repo_root; io = devnull)
end

# Spec a simular (ruta desde la raíz del repositorio)
spec   = joinpath(@__DIR__, "..", "config", "specs", "ejemplo_leo.json")

# Carpeta de resultados: con nothing, output/<nombre de la spec>/ (p. ej. output/ejemplo_leo/)
#output = joinpath(@__DIR__, "..", "output", "mi_prueba")
output = nothing

# main.jl solo se carga la primera vez en la sesión
if !isdefined(Main, :run_spec)
    include(joinpath(@__DIR__, "main.jl"))
end

run_spec(spec, output)
