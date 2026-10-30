module Parameters

# =============================================================================
# Módulo Parameters (parámetros maestros)
#
# Parámetros de configuración del programa: valores que fijamos los
# desarrolladores y que el usuario NO elige en el JSON. Se separan de
# Constants.jl, que contiene solo constantes físicas (μ, R, ω⊕...).
#
# Para cambiar el comportamiento del simulador en todos los escenarios, se
# edita aquí; para un escenario concreto, se usa la spec (config/specs/).
#
# Parámetros:
#   RANDOM_W_RANGE_DEG_S : rango [mín, máx] de |ω| (°/s) del modo de actitud "random"
#
# Dependencias: ninguna (se puede incluir antes que cualquier otro módulo).
# =============================================================================

export RANDOM_W_RANGE_DEG_S

# -----------------------------------------------------------------------------
# Actitud
# -----------------------------------------------------------------------------

# Rango de |ω| (°/s) del modo "random"; se sortea uniforme dentro de él.
# Acordado con el tutor (0.5–5 °/s), del orden de lo observado en basura espacial en LEO.
const RANDOM_W_RANGE_DEG_S = (0.5, 5.0)

end # module Parameters
