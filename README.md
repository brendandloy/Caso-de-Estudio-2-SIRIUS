# SIRIUS – Simulador de imágenes de objetos residentes en el espacio

**SIRIUS** (*Space Imaging for Relative Identification of Uncooperative targetS*) es una herramienta en Julia que simula las imágenes de objetos residentes en el espacio (RSOs) captadas por cámaras. Propaga observadores y objetivos, calcula la magnitud aparente de los objetivos y genera imágenes con campo estelar, PSF y ruido de sensor.

**Estado:** en desarrollo.

## Contexto

Este repositorio continúa SIRIUS dentro del **Caso de Estudio 2** del Máster Universitario en Sistemas Espaciales (MUSE, UPM). Parte del trabajo de casos de estudio anteriores:

- *Herramienta de simulación orbital para el estudio de basura espacial*: Ignacio Ricardo Ferrando Magdalena, CE2, enero de 2026. Arquitectura modular, propagación orbital y entorno.
- *Simulador fotométrico de basura espacial basado en el modelo de reflexión de Phong*: Alejandro Ortiz Tello, CE2, enero de 2026. Modelo fotométrico.
- *Simulación de cámara orbital para observación de RSOs*: Ignacio Ricardo Ferrando Magdalena, CE3, julio de 2026. Simulación de cámara.

El código heredado se subió en el commit `Inicialización Repositorio`.

Objetivos de esta iteración:

- **Definición de escenarios:** generador automático de configuraciones (módulo `Scenario`). En curso.
- **Dinámica de actitud de los objetos:** modos de actitud de los objetivos y de los satélites observadores, seguimiento de un RSO y montaje de la cámara (módulo `Attitude`). En curso.
- **Forma de los objetivos:** bloque `shape` con medidas concretas o sorteadas según la población real (módulo `Geometry`). Completado; paneles solares e importación de modelos STP/STL/OBJ, pendientes.
- **Mayor realismo en la magnitud relativa de los objetos de la imagen.** No iniciado.

El detalle de lo añadido está en [CHANGELOG.md](CHANGELOG.md).

## Requisitos e instalación

- Julia 1.13 (el `Manifest.toml` está generado con 1.13.1).
- Kernels SPICE y catálogo estelar en `data/`:
  - `de440s.bsp` (efemérides planetarias)
  - `naif0012.tls.txt` (segundos intercalares)
  - `I_239_vmag10_reduced.fits` (catálogo Hipparcos reducido a V ≤ 10)

Instalación de dependencias, desde la raíz del repositorio:

```julia
julia --project=.
] instantiate
```

Para comprobar la instalación, todas las pruebas de `test/` se ejecutan con `julia --project=. test/runtests.jl` (o `] test` en el REPL con el entorno activo).

El catálogo reducido se genera a partir del Hipparcos completo con `scripts/reduce_hipparcos_catalog.py` (Python + `astropy`).

## Uso rápido

El escenario se escribe a mano como una *spec* corta en `config/specs/`. `main.jl` genera primero el JSON completo en `config/generated/<nombre de la spec>.json` (lo sobrescribe si existe) y después lanza la simulación con él. La spec es obligatoria; la carpeta de resultados es opcional (por defecto, `output/<nombre de la spec>/`, p. ej. `output/ejemplo_leo/`, para no mezclar resultados de escenarios distintos). El JSON generado guarda las semillas usadas, así que una ejecución con elementos aleatorios se puede repetir.

Hay tres formas de ejecutarlo (desde la raíz del repositorio):

| Forma | Cómo | Velocidad |
|---|---|---|
| Terminal | `julia --project=. scripts/main.jl config/specs/ejemplo_leo.json [carpeta_resultados]` | Cada ejecución arranca un Julia nuevo y vuelve a cargar todos los paquetes |
| Lanzador de VS Code | Cambiar la variable `spec` (y, si se quiere, `output`) de `scripts/run_simulation_vscode.jl` y ejecutarlo con *Julia: Execute active File in REPL* (activa solo el entorno del repositorio) | La primera vez carga todo; las siguientes, en la misma sesión, van directas a la simulación |
| REPL a mano | `include("scripts/main.jl")` una vez y después `run_spec("config/specs/ejemplo_leo.json")` (o con una carpeta como segundo argumento) | Igual que el lanzador |

Para varias ejecuciones seguidas conviene el REPL (lanzador o `run_spec`); para una ejecución suelta, la terminal.

> **Aviso:** una sesión de REPL abierta sigue usando el código de `src/` que cargó al principio. Si se modifica algún archivo de `src/`, hay que cerrar el REPL (`exit()` o la papelera del terminal de Julia) y abrir uno nuevo; si no, se ejecuta el código antiguo.

La primera ejecución tras clonar el repositorio o cambiar `src/` tarda más porque Julia precompila (`julia --project=. -e "import Pkg; Pkg.precompile()"` lo hace por adelantado). Otros usos: simular directamente un JSON completo (lo hace `MonteCarlo.jl`) con `include("scripts/main.jl"); run_simulation("config/generated/ejemplo_leo.json", "output")`, y generar el JSON sin simular con `julia --project=. scripts/generate_config.jl <spec.json>`.

Ejemplo de *spec* (observador en LEO polar, objetivo en órbita kepleriana):

```json
{
  "scenario":   {"init_epoch": "2000-01-01T11:59:28", "Duration_hours": 0.2, "Data_time_step_seconds": 2},
  "objectives": {"deb":  {"orbit": {"type": "keplerian", "SMA": 8800, "ECC": 0.2, "INC": 90, "RAAN": 10, "AOP": 0, "TA": 357.5},
                          "attitude": {"type": "spin", "w_deg_s": 0},
                          "shape":    {"type": "cylinder", "radius_m": 1.5, "height_m": 8.0}}},
  "observers":  {"cam1": {"orbit": {"type": "LEO", "altitude_km": 620, "INC": 90, "RAAN": 0, "TA": 0},
                          "attitude": {"type": "pointing",
                                       "primary":   {"body_axis": [0, 0, 1],  "direction": [0, 1, 0]},
                                       "secondary": {"body_axis": [0, -1, 0], "direction": [0, 0, 1]}},
                          "camera_mount": {"boresight": [0, 0, 1], "up": [0, -1, 0]}}}
}
```

Hay más ejemplos en `config/specs/` (SSO + GEO, formaciones, CW, Walker, TLE, órbitas aleatorias y modos de actitud, formas y una cámara en `tracking` en `ejemplo_actitud.json`).

### Crear una spec nueva: `config/specs/plantilla.json`

La plantilla contiene **todas** las opciones que admite una spec, escritas con su valor por defecto: el bloque `scenario` completo (tiempos, *flags*, `camera_config`, `photometry_config`, `stations_llh`), el bloque `member_fields` (errores de navegación de las cámaras), la `seed` y un objetivo y un observador de ejemplo. Para un escenario nuevo:

1. Copiar `plantilla.json` en `config/specs/` con otro nombre.
2. Cambiar solo lo necesario. Por ejemplo, para guardar las imágenes: `"image_generation_flag": true` (conviene subir `Data_time_step_seconds`, porque se genera un PNG por paso y por cámara).
3. Las claves que no se cambien pueden borrarse: el generador usa el mismo valor por defecto.

Una clave mal escrita en la spec, en `scenario` o en `member_fields` (p. ej. `"image_generation_flg"`) da un error que lista las claves admitidas. Sin `seed`, cada ejecución sortea una nueva y la imprime.

## Estructura del repositorio

```
├── config/
│   ├── specs/         # Specs de escenario escritas a mano; plantilla.json lista todas las opciones
│   ├── generated/     # Configuraciones completas generadas a partir de las specs
│   ├── tle/           # Archivos TLE de ejemplo
│   ├── legacy/        # Configuraciones heredadas escritas a mano (formato anterior, sin "attitude" ni "shape")
│   └── *.json         # Configuración del Kalman y de Monte Carlo
├── data/              # Kernels SPICE y catálogo Hipparcos
├── output/            # Resultados de simulación (ignorado por git)
├── test/
│   ├── runtests.jl                 # Ejecuta todas las pruebas (cada archivo en su propio módulo)
│   ├── test_attitude.jl            # Pruebas del módulo Attitude (convenciones del T6 y modos)
│   ├── test_camera.jl              # Pruebas de CameraSim (intrínsecos, proyección pinhole, FOV del sensor, PSF, flujo)
│   ├── test_geometry.jl            # Pruebas de Geometry (mallas, L_c, ejes, bloque shape, sorteos)
│   ├── test_integrators.jl         # Pruebas de Integrators (Simpson exacto, orden 4, estado vectorial)
│   ├── test_photometry.jl          # Pruebas de Photometry (placa de una faceta, ley 1/d², sombra)
│   ├── test_propagator.jl          # Pruebas de Propagator (instantes de salida = malla 0:Δt:T)
│   ├── test_scenario.jl            # Pruebas del generador (reproducibilidad, elementos ↔ r-v, SSO)
│   ├── test_transformations.jl     # Pruebas de Transformations (Kepler, ERA, ECEF ↔ ECI, WGS84)
│   └── test_visibility.jl          # Pruebas de Visibility (eclipse, línea de visión, deslumbramiento)
├── scripts/
│   ├── main.jl                     # Lanzador de terminal: lee los argumentos y llama a run_spec (src/Simulation.jl)
│   ├── generate_config.jl          # Spec → JSON de escenario
│   ├── benchmark.jl                # Tiempo y memoria de la propagación y de la simulación completa de una spec
│   ├── Kalman_RSO.jl               # Filtro de Kalman sobre las medidas LOS simuladas
│   ├── MonteCarlo.jl               # Campañas Monte Carlo sobre la configuración
│   ├── run_*_vscode.jl             # Lanzadores para VS Code (en el de simulación se elige la spec)
│   └── reduce_hipparcos_catalog.py # Reducción del catálogo Hipparcos
└── src/
    ├── SIRIUS.jl           # Módulo raíz del paquete SIRIUS: incluye todos los submódulos
    ├── Types.jl            # Tipos de datos (State, TrajectoryData, AttitudeData…)
    ├── Constants.jl        # Constantes físicas (MU_EARTH, R_EARTH, J2_EARTH, OMEGA_EARTH, WGS84…)
    ├── Parameters.jl       # Parámetros maestros del programa (no los elige el usuario; ver abajo)
    ├── Integrators.jl      # Integradores de paso fijo (runge_kutta_4), genéricos para cualquier estado
    ├── JSON_utils.jl       # Lectores y validadores comunes de los bloques del JSON (claves, números, semillas)
    ├── Transformations.jl  # Rotaciones y conversiones de sistemas de referencia
    ├── Attitude/           # Actitud de objetivos y observadores y montaje de la cámara (ver abajo)
    │   ├── Attitude.jl       # Tipo base, registro de modos, lectura del bloque "attitude", bind_attitude
    │   ├── Quaternions.jl    # Cuaternión, Euler, DCM elementales, sorteo uniforme y propagación con RK4 (T6)
    │   ├── Frames.jl         # LVLH, direcciones de apuntado y método TRIAD
    │   ├── Mount.jl          # Bloque "camera_mount" (C_cb) y orientación de la cámara C_ci
    │   └── modes/            # spin, random, fixed, pointing, tracking
    ├── Visibility.jl       # Eclipse por la Tierra (umbra/penumbra), línea de visión y cono de deslumbramiento solar
    ├── Environment.jl      # Lectura de configuración, SPICE, dirección del Sol, error de navegación, exportación CSV
    ├── Propagator.jl       # Propagación orbital
    ├── Geometry/           # Forma de los objetivos (ver abajo)
    │   ├── Geometry.jl       # Bloque "shape", registro de formas, sorteo, ópticas y generate_debris
    │   ├── Meshes.jl         # Mallas, normales y áreas; envolvente convexa del fragmento
    │   ├── Population.jl     # Modelo de población (clases y tamaños de ESA) para los sorteos
    │   └── shapes/           # sphere, prism, plate, cylinder, cone, trunk, fragment
    ├── Photometry.jl       # Magnitud aparente (modelo de Phong por facetas)
    ├── Camera_sim.jl       # Generación de imagen con la orientación de la cámara: estrellas, PSF, fondo, ruido
    ├── Visualization.jl    # Visualización 3D de órbitas (GLMakie)
    ├── Simulation.jl       # Flujo completo: run_spec (spec → JSON → simulación) y run_simulation, por pasos
    └── Scenario/           # Generador de escenarios (ver abajo)
        ├── Scenario.jl       # Registros de plugins y autocarga
        ├── Builder.jl        # Valores por defecto, validación y escritura del JSON
        ├── Orbital_utils.jl  # Utilidades orbitales (SSO, LTAN, Kepler, r-v ↔ elementos)
        ├── Relative_utils.jl # Movimiento relativo (LVLH → elementos)
        ├── Orbits/           # keplerian, LEO, SSO, GEO, random_LEO
        ├── Patterns/         # single, formation, walker, tle
        └── Layouts/          # along_track, element_offsets, cw_state, cw_ellipse, pco
```

### Generador de escenarios (`Scenario`)

Convierte una *spec* corta en el JSON completo que lee `Environment.load_config`. Los valores que la spec no indica se toman de `DEFAULT_TEMPLATE` y `DEFAULT_MEMBER_FIELDS` (`Builder.jl`), los mismos que lista `config/specs/plantilla.json`. Se basa en tres tipos de *plugin*, cada uno en su carpeta:

| Tipo | Qué define | Tipos disponibles |
|---|---|---|
| `Orbits/` | Elementos keplerianos de un objeto | `keplerian`, `LEO`, `SSO`, `GEO`, `random_LEO` |
| `Patterns/` | Cómo se generan los miembros de un grupo | `single`, `formation`, `walker`, `tle` |
| `Layouts/` | Disposición de una formación respecto al *chief* | `along_track`, `element_offsets`, `cw_state`, `cw_ellipse`, `pco` |

Para añadir un tipo nuevo basta con crear un archivo en la carpeta correspondiente y registrarlo con `register_orbit!`, `register_pattern!` o `register_layout!`. No hay que modificar nada más.

Cada grupo de objetivos debe definir su actitud con un bloque `attitude` y su forma con un bloque `shape`; cada grupo de observadores, su actitud con `attitude` y, opcionalmente, el montaje de su cámara con `camera_mount`. Se aplican a todos los miembros del grupo. En el modo de actitud `random` y en las formas, cada miembro recibe su propia semilla a partir de la `seed` de la spec; lo que la forma no indica se sortea por miembro y se escribe ya resuelto (medidas concretas y `seed`) en el JSON generado.

### Actitud de objetivos y observadores (`Attitude`)

Cada objetivo y cada observador del JSON del simulador lleva un bloque `attitude` obligatorio con un modo (`type`). En los observadores describe la actitud del satélite, no la de la cámara. El módulo lo valida y calcula la DCM cuerpo → ECI en cada instante:

| Modo | Parámetros | Comportamiento |
|---|---|---|
| `spin` | `w_deg_s`, `w_dir` (ejes cuerpo), orientación inicial opcional | Giro con ω constante alrededor de un eje fijo en el cuerpo |
| `random` (solo objetivos) | `seed` (opcional; el generador pone una por miembro) | Orientación inicial uniforme, eje uniforme y \|ω\| uniforme en el rango del parámetro maestro `RANDOM_W_RANGE_DEG_S`; luego como `spin` |
| `fixed` | `frame` (`ECI` o `LVLH`), orientación | Orientación constante respecto a ECI o a la órbita |
| `pointing` | `primary` y `secondary`: `body_axis` + `direction` | Eje primario hacia `nadir`, `zenith`, `velocity`, `orbit_normal`, `sun` (con `-` para invertir) o un vector ECI; el secundario fija el giro (TRIAD) |
| `tracking` (solo observadores) | `objective`; `roll`: `body_axis` + `direction` | El boresight de la cámara apunta en cada instante al objetivo indicado; `roll` fija el giro del satélite alrededor de ese eje (TRIAD) |

La orientación se da como cuaternión `q` o como `euler_deg` + `sequence`. Para añadir un modo basta con crear un archivo en `src/Attitude/modes/` y registrarlo con `register_attitude!` (con `roles` si solo vale para objetivos u observadores). Detalle de cada modo: `?parse_spin`, `?parse_random`, `?parse_fixed`, `?parse_pointing`, `?parse_tracking`.

**Montaje de la cámara (`camera_mount`, observadores):** `boresight` (eje de visión) y `up` (dirección que sale arriba en la imagen), ambos en ejes cuerpo del satélite; `up` no puede ser paralelo a `boresight`. Sin bloque se usa `boresight = [0, 0, 1]`, `up = [0, -1, 0]` (ejes cámara = ejes cuerpo) y se avisa; el generador escribe el bloque en el JSON. Detalle: `?parse_camera_mount`, `?camera_dcm_history`. Los intrínsecos de la cámara no cambian: siguen saliendo de `fov_deg` y la resolución de `camera_config`.

### Forma de los objetivos (`Geometry`)

Cada objetivo lleva un bloque `shape` obligatorio con `type` y, opcionalmente, sus medidas en metros (todas o ninguna):

| `type` | Medidas | L_c = (X + Y + Z)/3 | Ejes cuerpo | Sin medidas, tamaño de |
|---|---|---|---|---|
| `sphere` | `radius_m` | 2r | Origen en el centro | objetos de misión |
| `prism` | `x_m`, `y_m`, `z_m` | (x + y + z)/3 | Origen en el centro, aristas según x, y, z; z longitudinal (sorteado: z ≥ x ≥ y) | satélites |
| `plate` | `x_m`, `y_m`, `thickness_m` | (x + y + t)/3 | Origen en el centro, z normal a la placa | objetos de misión |
| `cylinder` | `radius_m`, `height_m` | (4r + h)/3 | Origen en el centro, z eje de revolución | etapas de cohete |
| `cone` | `radius_m`, `height_m` | (4r + h)/3 | Origen en la base, z hacia el vértice | objetos de misión |
| `trunk` | `radius_m`, `top_radius_m`, `height_m` | (4R + h)/3 | Origen en la base mayor, z hacia la menor | objetos de misión |
| `fragment` | `lc_m` | `lc_m` | Origen en el centroide, ejes principales (x mayor) | fragmentos |

- **Medidas concretas:** se usan tal cual, sin mínimo. Si L_c < `SHAPE_LC_MIN_M` se avisa de que el objetivo probablemente no se verá en la imagen.
- **Sin medidas** (`{"type": "cylinder"}`): se sortean las proporciones de la forma y el tamaño de la clase indicada en la tabla.
- **`{"type": "random"}`:** se sortean la clase de objeto, la forma y las medidas según la población real (ver *Modelos físicos*); el JSON generado añade la `class` sorteada.
- **`seed`:** semilla de la malla del fragmento y de las propiedades ópticas; la pone el generador (una por miembro) y no se escribe en la spec.

Para añadir una forma basta con crear un archivo en `src/Geometry/shapes/` y registrarlo con `register_shape!` (contrato en la cabecera de `src/Geometry/Geometry.jl`). Detalle: `?parse_shape`, `?resolve_shape`, `?generate_debris`. Los paneles solares y la importación de modelos STP/STL/OBJ quedan para más adelante.

### Parámetros maestros (`Parameters`)

`src/Parameters.jl` reúne los parámetros de configuración del programa que fijamos los desarrolladores y que **no** se eligen en el JSON (a diferencia de `Constants.jl`, que solo tiene constantes físicas). Cambiarlos afecta a todos los escenarios. Contiene:

- `RANDOM_W_RANGE_DEG_S = (0.5, 5.0)`: rango de |ω| (°/s) del modo de actitud `random`.
- `ATTITUDE_MAX_STEP_DEG = 1.0`: ángulo máximo (°) que gira el cuerpo en un subpaso del RK4 de la actitud (ver Modelos físicos).
- `SHAPE_LC_MIN_M = 0.10`: L_c mínima (m) de los objetivos con medidas sorteadas.
- Modelo de población de las formas sorteadas: `SHAPE_CLASS_COUNTS`, `SHAPE_WEIGHTS_BY_CLASS`, `FRAGMENT_LC_BREAKS_M`, `FRAGMENT_LC_BETAS`, `PAYLOAD_AREA_BINS_M2`, `ROCKET_BODY_AREA_BINS_M2`, `MISSION_RELATED_LC_RANGE_M` y los rangos de proporciones de cada forma (`PRISM_RATIO_RANGE` … `TRUNK_RADIUS_RATIO_RANGE`). Cada uno indica su fuente en el archivo.
- Coeficientes ópticos de cada material: `OPTICS_BODY`, `OPTICS_TOP`, `OPTICS_PANEL` (medias y σ de absorción y reflexión especular, CE2 de Ortiz).
- Resolución de las mallas: `SPHERE_MESH_DIVISIONS = (20, 20)`, `CYLINDER_MESH_DIVISIONS = 20`, `CONE_MESH_DIVISIONS = 50` y `FRAGMENT_NUM_POINTS_RANGE = (10, 30)`.
- Escenarios: `RANDOM_LEO_ALTITUDE_RANGE_KM = (400, 2000)` (altura por defecto de `random_LEO`), `CW_MAX_ECC = 0.01` y `CW_MAX_RATIO = 1e-3` (avisos de validez de Clohessy–Wiltshire) y `TLE_MAX_DAYS = 3` (aviso por TLE antiguo).
- Propagación orbital: `PROPAGATOR_RELTOL = PROPAGATOR_ABSTOL = 1e-11` (tolerancias de `Vern9`).
- Cámara: `FOV_PREFILTER_MARGIN = 1.25` (margen del cono de prefiltro sobre la esquina del sensor).

Criterio: aquí van los valores que cambian el modelo o los resultados y que se podrían querer ajustar; las tolerancias que solo protegen un cálculo se quedan junto a su función.

## Modelos físicos

- **Propagación orbital:** gravedad central más armónicos zonales y teserales EGM2008 hasta grado y orden 4 (Montenbruck & Gill, 2000), integrada con `Vern9` (tolerancias `PROPAGATOR_RELTOL`/`PROPAGATOR_ABSTOL` = 1e-11). Unidades: km, km/s; sistema inercial ECI (J2000).
- **Actitud de los objetivos:** convenciones del tema T6 (*Attitude Kinematics*). Cuaternión escalar primero, `C(q) = (q_η² − q_vᵀq_v) I + 2 q_v q_vᵀ − 2 q_η [q_v×]` como DCM referencia → cuerpo; Euler con secuencia `abc` y `C = C_c(θ3) C_b(θ2) C_a(θ1)`. El simulador usa la DCM cuerpo → ECI (`AttitudeData.DCM_b2i`, la transpuesta).
  - `spin`/`random`: se propaga el cuaternión con la cinemática del T6, `q̇ = ½ Ω(ω) q`, integrada con Runge-Kutta 4 (`?runge_kutta_4`, `?quat_step`). Cada intervalo de salida Δt se divide en `n = ⌈|ω| Δt / θ_max⌉` subpasos, con θ_max = `ATTITUDE_MAX_STEP_DEG`, porque el error del RK4 depende del ángulo girado por paso y no del tiempo; con 1° el error de fase es ≈ 8e-13 rad por subpaso (≈ 3e-10 rad por vuelta). El cuaternión se renormaliza tras cada subpaso, porque el RK4 no conserva |q| = 1. Se usa RK4 y no la solución exacta para ω constante, `exp(½ Ω Δt)`, porque sirve igual cuando ω varíe (dinámica con inercia). La DCM cuerpo → ECI de cada instante es `C(q)ᵀ`. Sin inercia: el giro con ω constante solo es dinámicamente consistente si ê es un eje principal de inercia (hipótesis aceptada con el tutor, porque los objetos se ven como puntos).
  - Actitud inicial con ángulos de Euler: la DCM se pasa a cuaternión con la fórmula del T6 (diap. 36), `q_η = ½√(1 + tr C)`, `q_v = [C23 − C32, C31 − C13, C12 − C21]/(4 q_η)`. Como divide por q_η = cos(φ/2), un giro a menos de ≈ 0,1° de 180° (q_η < 1e-3) da error y hay que escribirlo como `q`, que se usa sin convertir (`?dcm_to_quat`).
  - `random`: orientación uniforme (cuaternión a partir de cuatro normales), eje uniforme en la esfera y |ω| uniforme en `RANDOM_W_RANGE_DEG_S` (0,5–5 °/s).
  - LVLH (igual que en `Relative_utils.jl`): x radial, z normal al plano orbital (r × v), y = z × x.
  - `pointing`: método TRIAD con r, v (ECI) y la dirección del Sol (geocéntrica, SPICE). Error si las dos direcciones son paralelas en algún instante.
  - `tracking`: TRIAD con el boresight de la cámara en ejes cuerpo hacia la línea de visión û = (r_obj − r_obs)/|r_obj − r_obs| y `roll.body_axis` hacia `roll.direction`. Supuestos: posiciones reales (no las de navegación `nav_`), sin tiempo de luz ni aberración (≈ 5″ en LEO, frente a ≈ 95″ por píxel con 27° y 1024 px); se apunta al objetivo aunque esté eclipsado u oculto.
- **Entorno:** dirección del Sol con SPICE (DE440s), eclipses por la Tierra, línea de visión y cono de deslumbramiento solar.
- **Geometría de los objetos:** cuerpos triangulados (prismas, esferas, placas, cilindros, conos, troncos y fragmentos) en ejes cuerpo, con medidas en m y mallas en km (áreas en km², las de la fotometría). Conos y troncos: solo superficie lateral exterior, sin bases (superficies abiertas). Fragmentos: envolvente convexa de 10–30 puntos aleatorios, en sus ejes principales y escalada a la longitud característica del NASA Standard Breakup Model, L_c = (X + Y + Z)/3. Sin paneles solares (pendiente).
  - **Ópticas:** coeficientes C_abs, C_spec, C_diff con normales truncadas alrededor de valores nominales (CE2 de Ortiz); el prisma pondera cuerpo y tapas ±z por área, la placa usa el material de panel. Todo lo aleatorio de un objetivo usa su `seed`.
  - **Modelo de población** (formas sorteadas), con datos de ESA's Annual Space Environment Report, issue 10.1 (2026): la clase de objeto se sortea en proporción a los objetos catalogados en LEO (tabla 3.1: 59,5 % satélites, 33,9 % fragmentos, 3,8 % etapas, 2,8 % objetos de misión con L_c ≥ 10 cm); los fragmentos siguen N(>L_c) ∝ L_c^−β por tramos (β = 2,17 de 1 mm a 1 cm y 1,36 de 1 a 10 cm, de MASTER 2026; β = 1,6 de 10 cm a 1 m, NASA SBM) y, si L_min < 10 cm, su peso crece con esa ley; el área media de satélites y etapas sale de los histogramas de las figuras 3.1b y 3.2b (lectura aproximada) y se convierte en tamaño con el teorema de Cauchy (área proyectada media = superficie/4, cuerpo convexo); los objetos de misión tienen L_c log-uniforme en 0,1–3 m (supuesto). Las proporciones de cada forma son uniformes en rangos aproximados de objetos reales. Supuestos: clase → forma (satélite → prisma, etapa → cilindro, fragmento → fragmento, objeto de misión → placa, cono o tronco) y clase independiente del tamaño por encima de 10 cm; el área de los satélites incluye paneles, así que el prisma equivalente sale mayor que el cuerpo real.
- **Fotometría:** modelo de Phong por facetas (componentes difusa y especular) para obtener la magnitud aparente V.
- **Orientación de la cámara:** ejes cámara dextrógiros (x_c derecha de la imagen, y_c abajo, z_c = boresight), coherentes con el T6 (triedros dextrógiros, det C = +1). `C_cb` sale de `camera_mount` (filas: x_c, y_c, z_c en ejes cuerpo) y `C_ci = C_cb · C_bi` en cada instante; su fila 3 es el boresight en ECI. Proyección pinhole (= gnomónica TAN, CE3 de I. Ferrando) con los intrínsecos de `camera_config` y la matriz K de la *Memoria Nagoya* (ec. 2.5): coordenadas normalizadas `x' = x_c/z_c`, `y' = y_c/z_c`, distorsión de Brown–Conrady (Brown, 1966 y 1971; el modelo de OpenCV y de `cameraIntrinsics` de MATLAB) `x'' = x' R + 2p1 x'y' + p2 (r² + 2x'²)`, `y'' = y' R + p1 (r² + 2y'²) + 2p2 x'y'`, con `R = 1 + k1 r² + k2 r⁴ + k3 r⁶` (`?distort`), y píxeles `u = f_x x'' + α y'' + u0`, `v = f_y y'' + v0` (`?radec_to_pixels!`). Las direcciones fuera del dominio de validez de la distorsión radial (donde `r R(r)` deja de crecer y la imagen se doblaría) no se proyectan. Claves: `resolution_x`, `resolution_y`, `fx_px`, `fy_px` (focales en px por radián; obligatorias), `u0_px`, `v0_px` (punto principal; `null` u omitidos = centro, (N + 1)/2) y `skew_px` (α, por defecto 0), y el bloque `distortion` con `k1`, `k2`, `k3`, `p1`, `p2` (cada uno por defecto 0; los mismos que `RadialDistortion` y `TangentialDistortion` de MATLAB). Píxeles desde 1 (el píxel k ocupa [k − ½, k + ½]), como en MATLAB; para intrínsecos de OpenCV, que cuenta desde 0, hay que sumar 1 a u0 y v0. El FOV ya no es una entrada: una fuente se ve si su píxel cae dentro del sensor (FOV rectangular, con las esquinas) y el FOV se calcula e imprime (`?sensor_fov_deg`); con el punto principal centrado, `FOV_x = 2 atan(N_x/2f_x)`. Los valores por defecto (`f_x = f_y = 512/tan(13,5°)`) dan el FOV de 27° del antiguo `fov_deg`. Antes de proyectar, el catálogo se filtra con un cono que contiene todo el sensor, con las esquinas pasadas a dirección deshaciendo K y la distorsión (`?prefilter_cone_deg`, `?undistort`). El FOV impreso es el nominal de K, sin distorsión. La PSF es una gaussiana con FWHM `fwhm_deg` pasada a píxeles con la escala del centro de la imagen en cada eje (`f_x · fwhm`, `f_y · fwhm`): circular con píxeles cuadrados y elíptica si `f_x ≠ f_y`. Con actitud constante en ECI el campo estelar se calcula una vez; si cambia, en cada imagen. La proyección anterior (RA/Dec con el norte arriba y el este a la derecha) daba la imagen en espejo (det = −1); ahora, con el norte arriba, el este sale a la izquierda, como en un sensor real, y los píxeles de `debris_positions.csv` cambian respecto a versiones anteriores.
- **Cámara:** catálogo Hipparcos dentro del campo de visión, PSF gaussiana, fondo, ruido Poisson, ruido de lectura y *bias*. Las medidas LOS se exportan con error. Por cada observador se escriben siempre `los_measurements.csv`, `state_measurements.csv` y `debris_positions.csv` (RSOs en el campo de visión con su posición en píxeles y su magnitud V en cada instante); los PNG solo con `image_generation_flag = true`.
- **Escenarios:**
  - Órbitas SSO con inclinación calculada por J2 y RAAN a partir de la LTAN.
  - Formaciones con condiciones iniciales de Clohessy-Wiltshire, válidas para un *chief* casi circular (e < 0.01) y separaciones ≪ SMA. Después se propagan con la dinámica completa.
  - Constelaciones Walker T/P/F (*delta* y *star*).
  - TLE propagados con SGP4 y convertidos de TEME a J2000 (IAU-76/FK5, sin EOP).

## Validación y resultados

- `test/test_attitude.jl` (109 comprobaciones superadas con Julia 1.13.1): matrices de Euler 3-2-1 y 3-1-3 y cuaternión frente a las expresiones del T6; conversión DCM → cuaternión con la fórmula del T6 (error claro en giros de 180°; en `spin`, un giro de 180° dado como `q` se usa tal cual); la propagación con RK4 coincide con la solución analítica `C_ib(0) R(ê, ωt)` durante un día con paso de 2 s a 4 °/s (diferencia < 1e-6; ≈ 1,5e-7 según el análisis de error) y en un solo intervalo de 120 s a 2 °/s, 240 subpasos (< 1e-9), el número de subpasos es ⌈|ω|Δt/θ_max⌉ y se mantiene |q| = 1 tras 10⁵ pasos; el modo `spin` reproduce el antiguo `propagate_attitude` (diferencia < 1e-7; ≈ 2e-8 esperado) y cumple `dC_ib/dt = C_ib [ω×]`; |ω| de `random` dentro del rango maestro y `random` rechazado en los observadores y `tracking` en los objetivos; media de las DCM aleatorias ≈ 0 (orientación uniforme); `fixed` en LVLH y `pointing` alinean los ejes pedidos; `camera_mount` da C_cb ortogonal con det = +1 (identidad por defecto); en `tracking` el boresight apunta al objetivo (con la cámara en +z y en +x cuerpo) y `roll` queda en el plano pedido; con el norte arriba, el este sale a la izquierda de la imagen. Se ejecuta con `julia --project=. test/test_attitude.jl`.
- `test/test_geometry.jl` (602 comprobaciones superadas con Julia 1.12.7): área de la esfera ≈ 4πr² (−1,2 % por ser un poliedro inscrito), del prisma, y lateral del cono y del tronco frente a la analítica (< 1 %); envolvente convexa y L_c del fragmento igual a la pedida, con los ejes principales ordenados; L_c de las fórmulas frente a la malla de cada forma (< 2 %) y orientación de los ejes de cilindro y cono; errores de validación del bloque `shape` y aviso por debajo de `SHAPE_LC_MIN_M`; sorteos reproducibles, proporciones dentro de los rangos maestros, L_c ≥ 0,1 m y proporción de fragmentos y ley de potencias en 4000 objetivos `random`; `generate_debris` sin normales NaN para todas las formas. La envolvente convexa propia se comparó con `scipy.spatial.ConvexHull` en 50 fragmentos (misma área, volumen y número de caras). Se ejecuta con `julia --project=. test/test_geometry.jl`.
- `test/test_propagator.jl` (4 comprobaciones superadas con Julia 1.13.1): la trayectoria propagada tiene exactamente los instantes de la malla 0:Δt:T que usa `scripts/main.jl` para la dirección del Sol, también cuando T no es múltiplo de Δt (900 s con paso de 120 s → 8 instantes). Antes de la corrección del 2026-11-06 fallaba (9 instantes frente a 8). Se ejecuta con `julia --project=. test/test_propagator.jl`.
- `test/test_visibility.jl` (19 comprobaciones superadas con Julia 1.13.1): casos geométricos calculados a mano. Eclipse: en el eje de la sombra, umbra; en el lado del Sol, iluminado; a 7000 km detrás de la Tierra, puntos a ±1 km de los radios analíticos de los conos de umbra (≈ 6345,9 km) y penumbra (≈ 6411,0 km) caen en la región correcta. LOS: bloqueada con la Tierra en medio o con el segmento 10 km bajo la superficie, visible 10 km por encima; un extremo 50 m bajo la superficie se proyecta y uno 200 km bajo ella no tiene LOS. Deslumbramiento: ángulo Sol-objetivo de 10° y 90° frente a un cono de 15°. Se ejecuta con `julia --project=. test/test_visibility.jl`.
- `test/test_transformations.jl` (15 comprobaciones superadas con Julia 1.13.1): Kepler → cartesianas en órbita circular (r y v sobre los ejes), en perigeo y apogeo (a(1 ∓ e) y velocidad en el perigeo) y en una órbita general (energía −μ/2a, |h| = √(μa(1−e²)), inclinación y RAAN a partir de h); ERA en J2000 (280,46061837504°) y avance diario; rotación ECEF → ECI; WGS84 en el ecuador y en el polo (b = a(1 − f)).
- `test/test_photometry.jl` (6 comprobaciones superadas con Julia 1.13.1): placa de una faceta de 1 m² de frente al Sol y al observador, F = Φ (C_diff + C_spec) A / d²; al duplicar la distancia la magnitud aumenta 5·log10(2) ≈ 1,505; la normal girada con la DCM da el mismo resultado; faceta de espaldas y objetivo en umbra sin magnitud.
- `test/test_camera.jl` (42 comprobaciones superadas con Julia 1.13.1): RA/Dec ↔ vector; intrínsecos (punto principal por defecto en el centro, `null` = centro, errores si falta una focal, una focal no es positiva o la resolución no es entera; FOV del sensor de 27° con `f = 512/tan(13,5°)` y vertical con 768 px; cono de prefiltro); proyección pinhole en una imagen de 1001 × 1001 (boresight en el píxel 501, una dirección con tan α = 100/f a 100 px a la derecha y hacia abajo, detrás de la cámara fuera de la imagen, punto principal desplazado, `f_y = 2 f_x` y skew, cámara girada); distorsión (sin coeficientes no cambia nada; radial con k1 = 0,1 en x' = 0,2 da 0,2008, y en píxeles 702 en lugar de 701; tangenciales p1 y p2 en (0,1, 0,1) calculados a mano; `undistort` recupera el punto; fuera del dominio de validez no hay imagen; con barril el cono de prefiltro crece; lectura del bloque y error con una clave desconocida); selección por sensor (una estrella cerca de la esquina, fuera del antiguo cono, entra; otra fuera del sensor no); núcleo de la PSF normalizado y a mitad de altura a FWHM/2, también elíptico; una fuente conserva su flujo; 5 mag = factor 100; cono del FOV.
- `test/test_scenario.jl` (11 comprobaciones superadas con Julia 1.13.1): la misma spec con la misma semilla da el mismo JSON y otra semilla, otra órbita y otra forma; elementos → r, v → elementos; el generador y `kepler_to_cartesian` dan el mismo estado; inclinación SSO a 700 km ≈ 98,19°.
- `test/test_integrators.jl` (6 comprobaciones superadas con Julia 1.13.1): un paso con ẋ = 4t³ es exacto (regla de Simpson); con ẋ = −x, al dividir el paso entre 2 el error en t = 1 baja entre 15 y 18 veces (orden 4); el oscilador armónico vuelve a [1, 0] tras un periodo con 1000 pasos (< 1e-9) y el estado sigue siendo un `SVector`. Se ejecuta con `julia --project=. test/test_integrators.jl`.
- `test/runtests.jl` ejecuta todas: 815 comprobaciones superadas con Julia 1.13.1, también con `] test`.
- Rendimiento (`scripts/benchmark.jl`, Julia 1.13.1, sin contar la compilación), 2026-11-06: `ejemplo_actitud.json` (0,2 h, paso de 120 s, 6 objetivos, 2 cámaras, con imágenes) tarda 0,75 s (502 MiB reservados, 32 % en el recolector de basura); la propagación supone 1,2 ms. Según el perfil (`Profile`), ≈ 75 % del tiempo es guardar los PNG (copia de la imagen, percentil para el contraste y codificación) y ≈ 7 % el campo estelar con actitud variable; la fotometría queda por debajo del 5 %. `ejemplo_leo.json` (0,2 h, paso de 2 s, sin imágenes) tarda 14 ms. Tras pasar la fotometría a `DebrisModel` y al bucle por facetas, las reservas no cambian (3,51 M, 501 MiB) y el tiempo de una ejecución a otra varía entre 0,75 y 1,26 s por la escritura de los PNG; la diferencia no se debe al cambio.
- Simulación completa de `config/specs/ejemplo_actitud.json` sin errores (≈ 10 s): seis objetivos con los cuatro modos de actitud y formas con medidas concretas, sorteadas y `random`.
- Orientación de la cámara, con `ejemplo_actitud.json` y propagación a dos cuerpos (sin `DifferentialEquations` ni PNG): `cam1` (apuntado fijo a ECI +y, norte arriba) da las 657 posiciones de RSOs y las 1818 estrellas de la versión anterior, con la misma fila y magnitud y la columna reflejada (x → 1025 − x, ±1 px por redondeo); `cam_track` (`tracking` a `sat_nadir`) ve a su objetivo en los 361 instantes, siempre en el píxel central (512–513, 512). Las 9 specs de `config/specs/` se generan sin errores.
- Resto de la herramienta: [PENDIENTE]

## Créditos

- Autor de esta iteración: Brendan de la Oliva Yanguas.
- Autores de iteraciones anteriores: Ignacio Ricardo Ferrando Magdalena y Alejandro Ortiz Tello.
- Tutor: Ángel Luis Porras Hermoso.
- Máster Universitario en Sistemas Espaciales (MUSE), ETSIAE – Universidad Politécnica de Madrid.
