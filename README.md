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
│   └── *.json         # Configuraciones heredadas (escritas a mano), Kalman y Monte Carlo
├── data/              # Kernels SPICE y catálogo Hipparcos
├── output/            # Resultados de simulación (ignorado por git)
├── test/
│   ├── test_attitude.jl            # Pruebas del módulo Attitude (convenciones del T6 y modos)
│   └── test_geometry.jl            # Pruebas de Geometry (mallas, L_c, ejes, bloque shape, sorteos)
├── scripts/
│   ├── main.jl                     # Spec → JSON completo → simulación (propagación → forma y actitud → fotometría → cámara → CSV)
│   ├── generate_config.jl          # Spec → JSON de escenario
│   ├── Kalman_RSO.jl               # Filtro de Kalman sobre las medidas LOS simuladas
│   ├── MonteCarlo.jl               # Campañas Monte Carlo sobre la configuración
│   ├── run_*_vscode.jl             # Lanzadores para VS Code (en el de simulación se elige la spec)
│   └── reduce_hipparcos_catalog.py # Reducción del catálogo Hipparcos
└── src/
    ├── Project.jl          # Módulo raíz: incluye todos los submódulos
    ├── Types.jl            # Tipos de datos (State, TrajectoryData, AttitudeData…)
    ├── Constants.jl        # Constantes físicas (μ, R, ω⊕, WGS84…)
    ├── Parameters.jl       # Parámetros maestros del programa (no los elige el usuario; ver abajo)
    ├── Transformations.jl  # Rotaciones y conversiones de sistemas de referencia
    ├── Attitude/           # Actitud de objetivos y observadores y montaje de la cámara (ver abajo)
    │   ├── Attitude.jl       # Tipo base, registro de modos, lectura del bloque "attitude", bind_attitude
    │   ├── Quaternions.jl    # Cuaternión, Euler, DCM elementales, sorteo uniforme y propagación (T6)
    │   ├── Frames.jl         # LVLH, direcciones de apuntado y método TRIAD
    │   ├── Mount.jl          # Bloque "camera_mount" (C_cb) y orientación de la cámara C_ci
    │   └── modes/            # spin, random, fixed, pointing, tracking
    ├── Environment.jl      # Lectura de configuración, SPICE, Sol, eclipses, LOS, exportación CSV
    ├── Propagator.jl       # Propagación orbital
    ├── Geometry/           # Forma de los objetivos (ver abajo)
    │   ├── Geometry.jl       # Bloque "shape", registro de formas, sorteo, ópticas y generate_debris
    │   ├── Meshes.jl         # Mallas, normales y áreas; envolvente convexa del fragmento
    │   ├── Population.jl     # Modelo de población (clases y tamaños de ESA) para los sorteos
    │   └── shapes/           # sphere, prism, plate, cylinder, cone, trunk, fragment
    ├── Photometry.jl       # Magnitud aparente (modelo de Phong por facetas)
    ├── Camera_sim.jl       # Generación de imagen con la orientación de la cámara: estrellas, PSF, fondo, ruido
    ├── Visualization.jl    # Visualización 3D de órbitas (GLMakie)
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
| `random` | `seed` (opcional; el generador pone una por miembro) | Orientación inicial uniforme, eje uniforme y \|ω\| uniforme en el rango del parámetro maestro `RANDOM_W_RANGE_DEG_S`; luego como `spin` |
| `fixed` | `frame` (`ECI` o `LVLH`), orientación | Orientación constante respecto a ECI o a la órbita |
| `pointing` | `primary` y `secondary`: `body_axis` + `direction` | Eje primario hacia `nadir`, `zenith`, `velocity`, `orbit_normal`, `sun` (con `-` para invertir) o un vector ECI; el secundario fija el giro (TRIAD) |
| `tracking` (solo observadores) | `objective`; `roll`: `body_axis` + `direction` | El boresight de la cámara apunta en cada instante al objetivo indicado; `roll` fija el giro del satélite alrededor de ese eje (TRIAD) |

La orientación se da como cuaternión `q` o como `euler_deg` + `sequence`. Para añadir un modo basta con crear un archivo en `src/Attitude/modes/` y registrarlo con `register_attitude!`. Detalle de cada modo: `?parse_spin`, `?parse_random`, `?parse_fixed`, `?parse_pointing`, `?parse_tracking`.

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
- `SHAPE_LC_MIN_M = 0.10`: L_c mínima (m) de los objetivos con medidas sorteadas.
- Modelo de población de las formas sorteadas: `SHAPE_CLASS_COUNTS`, `SHAPE_WEIGHTS_BY_CLASS`, `FRAGMENT_LC_BREAKS_M`, `FRAGMENT_LC_BETAS`, `PAYLOAD_AREA_BINS_M2`, `ROCKET_BODY_AREA_BINS_M2`, `MISSION_RELATED_LC_RANGE_M` y los rangos de proporciones de cada forma (`PRISM_RATIO_RANGE` … `TRUNK_RADIUS_RATIO_RANGE`). Cada uno indica su fuente en el archivo.

## Modelos físicos

- **Propagación orbital:** gravedad central más armónicos zonales y teserales EGM2008 hasta grado y orden 4 (Montenbruck & Gill, 2000), integrada con `Vern9` (tolerancias 1e-11). Unidades: km, km/s; sistema inercial ECI (J2000).
- **Actitud de los objetivos:** convenciones del tema T6 (*Attitude Kinematics*). Cuaternión escalar primero, `C(q) = (q_η² − q_vᵀq_v) I + 2 q_v q_vᵀ − 2 q_η [q_v×]` como DCM referencia → cuerpo; Euler con secuencia `abc` y `C = C_c(θ3) C_b(θ2) C_a(θ1)`. El simulador usa la DCM cuerpo → ECI (`AttitudeData.DCM_b2i`, la transpuesta).
  - `spin`/`random`: se propaga el cuaternión con la cinemática del T6, `q̇ = ½ Ω(ω) q`, paso a paso: `q_{k+1} = exp(½ Ω Δt) q_k = [cos(|ω|Δt/2) I₄ + sin(|ω|Δt/2)/|ω| · Ω] q_k` (exacta para ω constante, porque `Ω² = −|ω|² I₄`), renormalizando en cada paso para eliminar el redondeo acumulado. La DCM cuerpo → ECI de cada instante es `C(q)ᵀ`. Sin inercia: el giro con ω constante solo es dinámicamente consistente si ê es un eje principal de inercia (hipótesis aceptada con el tutor, porque los objetos se ven como puntos).
  - `random`: orientación uniforme (cuaternión a partir de cuatro normales), eje uniforme en la esfera y |ω| uniforme en `RANDOM_W_RANGE_DEG_S` (0,5–5 °/s).
  - LVLH (igual que en `Relative_utils.jl`): x radial, z normal al plano orbital (r × v), y = z × x.
  - `pointing`: método TRIAD con r, v (ECI) y la dirección del Sol (geocéntrica, SPICE). Error si las dos direcciones son paralelas en algún instante.
  - `tracking`: TRIAD con el boresight de la cámara en ejes cuerpo hacia la línea de visión û = (r_obj − r_obs)/|r_obj − r_obs| y `roll.body_axis` hacia `roll.direction`. Supuestos: posiciones reales (no las de navegación `nav_`), sin tiempo de luz ni aberración (≈ 5″ en LEO, frente a ≈ 95″ por píxel con 27° y 1024 px); se apunta al objetivo aunque esté eclipsado u oculto.
- **Entorno:** dirección del Sol con SPICE (DE440s), eclipses por la Tierra, línea de visión y cono de deslumbramiento solar.
- **Geometría de los objetos:** cuerpos triangulados (prismas, esferas, placas, cilindros, conos, troncos y fragmentos) en ejes cuerpo, con medidas en m y mallas en km (áreas en km², las de la fotometría). Conos y troncos: solo superficie lateral exterior, sin bases (superficies abiertas). Fragmentos: envolvente convexa de 10–30 puntos aleatorios, en sus ejes principales y escalada a la longitud característica del NASA Standard Breakup Model, L_c = (X + Y + Z)/3. Sin paneles solares (pendiente).
  - **Ópticas:** coeficientes C_abs, C_spec, C_diff con normales truncadas alrededor de valores nominales (CE2 de Ortiz); el prisma pondera cuerpo y tapas ±z por área, la placa usa el material de panel. Todo lo aleatorio de un objetivo usa su `seed`.
  - **Modelo de población** (formas sorteadas), con datos de ESA's Annual Space Environment Report, issue 10.1 (2026): la clase de objeto se sortea en proporción a los objetos catalogados en LEO (tabla 3.1: 59,5 % satélites, 33,9 % fragmentos, 3,8 % etapas, 2,8 % objetos de misión con L_c ≥ 10 cm); los fragmentos siguen N(>L_c) ∝ L_c^−β por tramos (β = 2,17 de 1 mm a 1 cm y 1,36 de 1 a 10 cm, de MASTER 2026; β = 1,6 de 10 cm a 1 m, NASA SBM) y, si L_min < 10 cm, su peso crece con esa ley; el área media de satélites y etapas sale de los histogramas de las figuras 3.1b y 3.2b (lectura aproximada) y se convierte en tamaño con el teorema de Cauchy (área proyectada media = superficie/4, cuerpo convexo); los objetos de misión tienen L_c log-uniforme en 0,1–3 m (supuesto). Las proporciones de cada forma son uniformes en rangos aproximados de objetos reales. Supuestos: clase → forma (satélite → prisma, etapa → cilindro, fragmento → fragmento, objeto de misión → placa, cono o tronco) y clase independiente del tamaño por encima de 10 cm; el área de los satélites incluye paneles, así que el prisma equivalente sale mayor que el cuerpo real.
- **Fotometría:** modelo de Phong por facetas (componentes difusa y especular) para obtener la magnitud aparente V.
- **Orientación de la cámara:** ejes cámara dextrógiros (x_c derecha de la imagen, y_c abajo, z_c = boresight), coherentes con el T6 (triedros dextrógiros, det C = +1). `C_cb` sale de `camera_mount` (filas: x_c, y_c, z_c en ejes cuerpo) y `C_ci = C_cb · C_bi` en cada instante; su fila 3 es el boresight en ECI. Proyección pinhole (= gnomónica TAN): `x = c_x + f x_c/z_c`, `y = c_y + f y_c/z_c`, con `f = ½ mín(N_x, N_y)/tan(FOV/2)` como en el CE3; FOV = cono de semiángulo FOV/2 alrededor del boresight. Con actitud constante en ECI el campo estelar se calcula una vez; si cambia, en cada imagen. La proyección anterior (RA/Dec con el norte arriba y el este a la derecha) daba la imagen en espejo (det = −1); ahora, con el norte arriba, el este sale a la izquierda, como en un sensor real, y los píxeles de `debris_positions.csv` cambian respecto a versiones anteriores.
- **Cámara:** catálogo Hipparcos dentro del campo de visión, PSF gaussiana, fondo, ruido Poisson, ruido de lectura y *bias*. Las medidas LOS se exportan con error. Por cada observador se escriben siempre `los_measurements.csv`, `state_measurements.csv` y `debris_positions.csv` (RSOs en el campo de visión con su posición en píxeles y su magnitud V en cada instante); los PNG solo con `image_generation_flag = true`.
- **Escenarios:**
  - Órbitas SSO con inclinación calculada por J2 y RAAN a partir de la LTAN.
  - Formaciones con condiciones iniciales de Clohessy-Wiltshire, válidas para un *chief* casi circular (e < 0.01) y separaciones ≪ SMA. Después se propagan con la dinámica completa.
  - Constelaciones Walker T/P/F (*delta* y *star*).
  - TLE propagados con SGP4 y convertidos de TEME a J2000 (IAU-76/FK5, sin EOP).

## Validación y resultados

- `test/test_attitude.jl` (100 comprobaciones superadas con Julia 1.13.3): matrices de Euler 3-2-1 y 3-1-3 y cuaternión frente a las expresiones del T6; conversión DCM ↔ cuaternión (incluidos giros de 180°); la propagación paso a paso coincide con la solución analítica `C_ib(0) R(ê, ωt)` durante un día con paso de 2 s (diferencia < 1e-9; ≈ 5e-13 en la práctica) y mantiene |q| = 1 tras 10⁵ pasos; el modo `spin` reproduce el antiguo `propagate_attitude` (diferencia < 1e-10) y cumple `dC_ib/dt = C_ib [ω×]`; |ω| de `random` dentro del rango maestro; media de las DCM aleatorias ≈ 0 (orientación uniforme); `fixed` en LVLH y `pointing` alinean los ejes pedidos; `camera_mount` da C_cb ortogonal con det = +1 (identidad por defecto); en `tracking` el boresight apunta al objetivo (con la cámara en +z y en +x cuerpo) y `roll` queda en el plano pedido; con el norte arriba, el este sale a la izquierda de la imagen. Se ejecuta con `julia --project=. test/test_attitude.jl`.
- `test/test_geometry.jl` (602 comprobaciones superadas con Julia 1.12.7): área de la esfera ≈ 4πr² (−1,2 % por ser un poliedro inscrito), del prisma, y lateral del cono y del tronco frente a la analítica (< 1 %); envolvente convexa y L_c del fragmento igual a la pedida, con los ejes principales ordenados; L_c de las fórmulas frente a la malla de cada forma (< 2 %) y orientación de los ejes de cilindro y cono; errores de validación del bloque `shape` y aviso por debajo de `SHAPE_LC_MIN_M`; sorteos reproducibles, proporciones dentro de los rangos maestros, L_c ≥ 0,1 m y proporción de fragmentos y ley de potencias en 4000 objetivos `random`; `generate_debris` sin normales NaN para todas las formas. La envolvente convexa propia se comparó con `scipy.spatial.ConvexHull` en 50 fragmentos (misma área, volumen y número de caras). Se ejecuta con `julia --project=. test/test_geometry.jl`.
- Simulación completa de `config/specs/ejemplo_actitud.json` sin errores (≈ 10 s): seis objetivos con los cuatro modos de actitud y formas con medidas concretas, sorteadas y `random`.
- Orientación de la cámara, con `ejemplo_actitud.json` y propagación a dos cuerpos (sin `DifferentialEquations` ni PNG): `cam1` (apuntado fijo a ECI +y, norte arriba) da las 657 posiciones de RSOs y las 1818 estrellas de la versión anterior, con la misma fila y magnitud y la columna reflejada (x → 1025 − x, ±1 px por redondeo); `cam_track` (`tracking` a `sat_nadir`) ve a su objetivo en los 361 instantes, siempre en el píxel central (512–513, 512). Las 9 specs de `config/specs/` se generan sin errores.
- Resto de la herramienta: [PENDIENTE]

## Créditos

- Autor de esta iteración: Brendan de la Oliva Yanguas.
- Autores de iteraciones anteriores: Ignacio Ricardo Ferrando Magdalena y Alejandro Ortiz Tello.
- Tutor: Ángel Luis Porras Hermoso.
- Máster Universitario en Sistemas Espaciales (MUSE), ETSIAE – Universidad Politécnica de Madrid.
