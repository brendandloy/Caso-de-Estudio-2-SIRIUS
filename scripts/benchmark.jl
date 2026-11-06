# =============================================================================
# benchmark.jl
#
# Mide el tiempo y la memoria reservada de la propagación y de la simulación
# completa de una spec, para comparar el rendimiento antes y después de un cambio.
#
# Uso (desde la raíz del repositorio):
#   julia --project=. scripts/benchmark.jl [<spec.json>]
#   (sin argumento usa config/specs/ejemplo_actitud.json)
#
# Cada medida se repite dos veces y solo se imprime la segunda: la primera llamada
# incluye la compilación de Julia, que no depende del código del proyecto. La
# simulación escribe en una carpeta temporal, para no tocar output/.
#
# Lo que imprime @time: tiempo, número de reservas de memoria ("allocations"),
# memoria total reservada y % del tiempo en el recolector de basura (gc time).
#
# Funciones: ninguna (llama a run_simulation y a propagate).
# =============================================================================

include(joinpath(@__DIR__, "main.jl"))      # carga SIRIUS; run_simulation

using SIRIUS.Scenario: generate_config_file
using SIRIUS.Environment: load_config
using SIRIUS.Propagator: propagate

# Spec a medir: la del argumento o la de ejemplo con actitud variable (tracking)
spec_path = length(ARGS) >= 1 ? abspath(ARGS[1]) :
            normpath(joinpath(@__DIR__, "..", "config", "specs", "ejemplo_actitud.json"))
spec_name = splitext(basename(spec_path))[1]

# JSON completo generado a partir de la spec (con semilla en la spec sale siempre igual)
cfg_path = normpath(joinpath(@__DIR__, "..", "config", "generated", spec_name * ".json"))
generate_config_file(spec_path, cfg_path)
env = load_config(cfg_path)
duration_s = env.Duration_hours * 3600                    # duración de la simulación (s)
step_s = env.Data_time_step_seconds                       # paso de salida (s)

# --- Propagación de todas las órbitas del escenario ---
propagate(env.objects, env.init_epoch, duration_s, step_s)          # primera vez: compila
println("\n=== propagate ($spec_name) ===")
@time propagate(env.objects, env.init_epoch, duration_s, step_s)

# --- Simulación completa (propagación, fotometría, cámara, CSV) ---
output_dir = mktempdir()                                  # carpeta temporal, se borra al salir
run_simulation(cfg_path, output_dir)                      # primera vez: compila
println("\n=== run_simulation ($spec_name) ===")
@time run_simulation(cfg_path, output_dir)
