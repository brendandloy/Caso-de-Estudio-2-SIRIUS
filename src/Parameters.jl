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
#   ATTITUDE_MAX_STEP_DEG       : ángulo máximo girado en un subpaso del RK4 de la actitud (°)
#   SHAPE_LC_MIN_M              : longitud característica mínima de los objetivos sorteados (m)
#   SHAPE_CLASS_COUNTS          : número de objetos catalogados por clase (ESA, LEO, > 10 cm)
#   SHAPE_WEIGHTS_BY_CLASS      : probabilidad de cada forma dentro de cada clase
#   FRAGMENT_LC_BREAKS_M,
#   FRAGMENT_LC_BETAS           : ley de potencias por tramos del tamaño de los fragmentos
#   PAYLOAD_AREA_BINS_M2,
#   ROCKET_BODY_AREA_BINS_M2    : histogramas de área media de satélites y etapas (ESA, LEO)
#   MISSION_RELATED_LC_RANGE_M  : rango de L_c de los objetos de misión (log-uniforme)
#   PRISM_RATIO_RANGE ... TRUNK_RADIUS_RATIO_RANGE : proporciones de cada forma
#   OPTICS_BODY, OPTICS_TOP, OPTICS_PANEL : distribución de los coeficientes ópticos por material
#   SPHERE_MESH_DIVISIONS ... FRAGMENT_NUM_POINTS_RANGE : resolución de las mallas
#   RANDOM_LEO_ALTITUDE_RANGE_KM : rango de altura por defecto de la órbita "random_LEO"
#   CW_MAX_ECC, CW_MAX_RATIO    : límites de validez de Clohessy-Wiltshire (avisos)
#   TLE_MAX_DAYS                : antigüedad máxima de un TLE respecto a init_epoch (aviso)
#   PROPAGATOR_RELTOL, PROPAGATOR_ABSTOL : tolerancias del integrador orbital (Vern9)
#   FOV_PREFILTER_MARGIN        : margen del cono de prefiltro de la cámara sobre la esquina del sensor
#
# Criterio: aquí van los valores que cambian el modelo o los resultados y que se
# podrían querer ajustar. Las tolerancias que solo protegen un cálculo (evitar
# dividir por ≈ 0, aceptar el redondeo de un dato) se quedan junto a su función.
#
# Fuente principal de la forma de los objetivos: ESA's Annual Space Environment
# Report, GEN-DB-LOG-00288-OPS-SD, Issue 10.1 (8 sep 2026), sección 2.4 (MASTER,
# población de referencia 01/02/2026), tabla 3.1 y figuras 3.1b y 3.2b.
#
# Dependencias: ninguna (se puede incluir antes que cualquier otro módulo).
# =============================================================================

export RANDOM_W_RANGE_DEG_S, ATTITUDE_MAX_STEP_DEG
export SHAPE_LC_MIN_M, SHAPE_CLASS_COUNTS, SHAPE_WEIGHTS_BY_CLASS
export FRAGMENT_LC_BREAKS_M, FRAGMENT_LC_BETAS
export PAYLOAD_AREA_BINS_M2, ROCKET_BODY_AREA_BINS_M2, MISSION_RELATED_LC_RANGE_M
export PRISM_RATIO_RANGE, PLATE_ASPECT_RANGE, PLATE_THICKNESS_RATIO_RANGE
export CYLINDER_HR_RANGE, CONE_HR_RANGE, TRUNK_HR_RANGE, TRUNK_RADIUS_RATIO_RANGE
export OPTICS_BODY, OPTICS_TOP, OPTICS_PANEL
export SPHERE_MESH_DIVISIONS, CYLINDER_MESH_DIVISIONS, CONE_MESH_DIVISIONS, FRAGMENT_NUM_POINTS_RANGE
export RANDOM_LEO_ALTITUDE_RANGE_KM, CW_MAX_ECC, CW_MAX_RATIO, TLE_MAX_DAYS
export PROPAGATOR_RELTOL, PROPAGATOR_ABSTOL
export FOV_PREFILTER_MARGIN

# -----------------------------------------------------------------------------
# Actitud
# -----------------------------------------------------------------------------

# Rango de |ω| (°/s) del modo "random"; se sortea uniforme dentro de él.
# Acordado con el tutor (0.5–5 °/s), del orden de lo observado en basura espacial en LEO.
const RANDOM_W_RANGE_DEG_S = (0.5, 5.0)

# Ángulo máximo (°) que gira el cuerpo en un subpaso del RK4 de la cinemática del cuaternión.
# El paso de salida (Data_time_step_seconds) se divide en n = ⌈|ω|Δt / θ_max⌉ subpasos, así la
# precisión no depende de Δt ni de |ω|. Con 1° el error de fase es ≈ 8e-13 rad por subpaso
# (≈ 3e-10 rad por vuelta completa).
const ATTITUDE_MAX_STEP_DEG = 1.0

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

# -----------------------------------------------------------------------------
# Propiedades ópticas y mallas de los objetivos
# -----------------------------------------------------------------------------

# (μ_abs, σ_abs, μ_spec_aux, σ_spec_aux) de cada material (valores heredados del CE2 de Ortiz)
const OPTICS_BODY  = (0.62, 0.12, 0.7, 0.1)        # cuerpo (MLI, aluminio...)
const OPTICS_TOP   = (0.55, 0.2, 0.7, 0.1)         # tapas ±z del prisma
const OPTICS_PANEL = (0.75, 0.1, 0.8, 0.05)        # paneles y placas

# Resolución de las mallas: más divisiones aproximan mejor el área y la normal de las
# superficies curvas (la esfera de 20 × 20 tiene un área un 1,2 % menor que 4πr²)
const SPHERE_MESH_DIVISIONS = (20, 20)             # (longitud nθ, latitud nφ con los polos)
const CYLINDER_MESH_DIVISIONS = 20                 # divisiones del contorno del cilindro
const CONE_MESH_DIVISIONS = 50                     # divisiones del contorno de conos y troncos
const FRAGMENT_NUM_POINTS_RANGE = (10, 30)         # puntos de la envolvente convexa de un fragmento

# -----------------------------------------------------------------------------
# Escenarios
# -----------------------------------------------------------------------------

# Rango de altura (km) por defecto de "random_LEO" si la spec no lo da (región LEO)
const RANDOM_LEO_ALTITUDE_RANGE_KM = (400.0, 2000.0)

# Clohessy-Wiltshire: excentricidad máxima del chief y cociente distancia relativa / SMA
# por encima de los cuales la aproximación pierde precisión (se avisa)
const CW_MAX_ECC = 0.01
const CW_MAX_RATIO = 1e-3

# Distancia máxima (días) entre la época de un TLE e init_epoch antes de avisar
const TLE_MAX_DAYS = 3.0

# -----------------------------------------------------------------------------
# Propagación orbital
# -----------------------------------------------------------------------------

# Tolerancias relativa y absoluta del integrador adaptativo Vern9 (Propagator.jl)
const PROPAGATOR_RELTOL = 1e-11
const PROPAGATOR_ABSTOL = 1e-11

# -----------------------------------------------------------------------------
# Cámara
# -----------------------------------------------------------------------------

# El catálogo se filtra primero con un cono alrededor del boresight cuyo semiángulo es el de
# la esquina del sensor más alejada por este factor; el margen cubre los puntos que la
# distorsión de la lente lleva al sensor desde fuera (distorsión de barril, k1 < 0)
const FOV_PREFILTER_MARGIN = 1.25

end # module Parameters
