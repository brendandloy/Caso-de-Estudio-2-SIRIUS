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
#   RANDOM_W_RANGE_DEG_S        : rango [mín, máx] de |ω| (°/s) del modo de actitud "random"
#   SHAPE_LC_MIN_M              : longitud característica mínima de los objetivos sorteados (m)
#   SHAPE_CLASS_COUNTS          : número de objetos catalogados por clase (ESA, LEO, > 10 cm)
#   SHAPE_WEIGHTS_BY_CLASS      : probabilidad de cada forma dentro de cada clase
#   FRAGMENT_LC_BREAKS_M,
#   FRAGMENT_LC_BETAS           : ley de potencias por tramos del tamaño de los fragmentos
#   PAYLOAD_AREA_BINS_M2,
#   ROCKET_BODY_AREA_BINS_M2    : histogramas de área media de satélites y etapas (ESA, LEO)
#   MISSION_RELATED_LC_RANGE_M  : rango de L_c de los objetos de misión (log-uniforme)
#   PRISM_RATIO_RANGE ... TRUNK_RADIUS_RATIO_RANGE : proporciones de cada forma
#
# Fuente principal de la forma de los objetivos: ESA's Annual Space Environment
# Report, GEN-DB-LOG-00288-OPS-SD, Issue 10.1 (8 sep 2026), sección 2.4 (MASTER,
# población de referencia 01/02/2026), tabla 3.1 y figuras 3.1b y 3.2b.
#
# Dependencias: ninguna (se puede incluir antes que cualquier otro módulo).
# =============================================================================

export RANDOM_W_RANGE_DEG_S
export SHAPE_LC_MIN_M, SHAPE_CLASS_COUNTS, SHAPE_WEIGHTS_BY_CLASS
export FRAGMENT_LC_BREAKS_M, FRAGMENT_LC_BETAS
export PAYLOAD_AREA_BINS_M2, ROCKET_BODY_AREA_BINS_M2, MISSION_RELATED_LC_RANGE_M
export PRISM_RATIO_RANGE, PLATE_ASPECT_RANGE, PLATE_THICKNESS_RATIO_RANGE
export CYLINDER_HR_RANGE, CONE_HR_RANGE, TRUNK_HR_RANGE, TRUNK_RADIUS_RATIO_RANGE

# -----------------------------------------------------------------------------
# Actitud
# -----------------------------------------------------------------------------

# Rango de |ω| (°/s) del modo "random"; se sortea uniforme dentro de él.
# Acordado con el tutor (0.5–5 °/s), del orden de lo observado en basura espacial en LEO.
const RANDOM_W_RANGE_DEG_S = (0.5, 5.0)

# -----------------------------------------------------------------------------
# Forma de los objetivos: tamaño mínimo y clases de objeto
# -----------------------------------------------------------------------------

# L_c mínima (m) de los objetivos con medidas sorteadas. 10 cm es el límite del catálogo
# en LEO; debajo casi todo son fragmentos y una cámara óptica difícilmente los ve.
const SHAPE_LC_MIN_M = 0.10

# Objetos catalogados en LEO a finales de 2025 por clase (ESA, tabla 3.1):
#   fragment        = PF + PD + RF + RD + UI (fragmentos, restos y no identificados)
#   payload         = PL (satélites, activos o no)
#   rocket_body     = RB (etapas de cohete)
#   mission_related = PM + RM (tapas, cofias, adaptadores, motores...)
const SHAPE_CLASS_COUNTS = (fragment = 8483, payload = 14906, rocket_body = 957, mission_related = 714)

# Probabilidad de cada forma dentro de su clase (supuesto nuestro; la ESA da clases, no formas).
# mission_related: PM (209, tapas) -> placa; RM (505) repartido entre placa, cono y tronco.
const SHAPE_WEIGHTS_BY_CLASS = (
    fragment        = (fragment = 1.0,),
    payload         = (prism = 1.0,),
    rocket_body     = (cylinder = 1.0,),
    mission_related = (plate = 0.53, cone = 0.235, trunk = 0.235),
)

# -----------------------------------------------------------------------------
# Distribuciones de tamaño de cada clase
# -----------------------------------------------------------------------------

# Fragmentos: N(>L_c) ∝ L_c^-β por tramos. Límites de los tramos (m):
const FRAGMENT_LC_BREAKS_M = (0.001, 0.01, 0.10, 1.0)

# β de cada tramo:
#   1 mm–1 cm : log10(231.6e6 / 1.568e6) = 2.17  (MASTER 2026, ESA §2.4)
#   1–10 cm   : log10(1.568e6 / 68 450)  = 1.36  (MASTER 2026, ESA §2.4)
#   10 cm–1 m : 1.6, NASA Standard Breakup Model para explosiones (Johnson et al., 2001)
const FRAGMENT_LC_BETAS = (2.17, 1.36, 1.6)

# Satélites: histograma del área media (m²) en LEO, (área mín, área máx, número de objetos).
# Lectura aproximada de la figura 3.1b del informe de ESA, agrupada en intervalos.
const PAYLOAD_AREA_BINS_M2 = (
    (0.001, 0.01, 40), (0.01, 0.1, 890), (0.1, 1.0, 1700),
    (1.0, 10.0, 3150), (10.0, 20.0, 3000), (20.0, 50.0, 6000),
)

# Etapas de cohete: igual, a partir de la figura 3.2b del informe de ESA.
const ROCKET_BODY_AREA_BINS_M2 = (
    (0.5, 1.0, 110), (1.0, 10.0, 590), (10.0, 20.0, 290), (20.0, 60.0, 460),
)

# Objetos de misión: sin datos de tamaño en el informe; L_c log-uniforme (supuesto nuestro).
const MISSION_RELATED_LC_RANGE_M = (0.10, 3.0)

# -----------------------------------------------------------------------------
# Proporciones de cada forma (se sortean uniformes dentro del rango)
# -----------------------------------------------------------------------------

# Prisma: las otras dos aristas entre 1 y 3 veces la menor (CubeSats 1U–6U y plataformas)
const PRISM_RATIO_RANGE = (1.0, 3.0)

# Placa: lado mayor / lado menor, y espesor / lado mayor (tapas, paneles)
const PLATE_ASPECT_RANGE = (1.0, 4.0)
const PLATE_THICKNESS_RATIO_RANGE = (0.005, 0.05)

# Cilindro: altura / radio (etapas en órbita: L/D ≈ 1.5–4, p. ej. SL-8, SL-16, Falcon 9 S2)
const CYLINDER_HR_RANGE = (3.0, 8.0)

# Cono: altura / radio de la base (toberas, adaptadores)
const CONE_HR_RANGE = (1.0, 4.0)

# Tronco de cono: altura / radio mayor, y radio mayor / radio menor (adaptadores de carga)
const TRUNK_HR_RANGE = (0.5, 2.0)
const TRUNK_RADIUS_RATIO_RANGE = (1.2, 3.0)

end # module Parameters
