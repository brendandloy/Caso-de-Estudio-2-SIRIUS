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
- **Dinámica de actitud de los objetos:** modos de actitud de los objetivos (módulo `Attitude`). En curso.
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
                          "attitude": {"type": "spin", "w_deg_s": 0}}},
  "observers":  {"cam1": {"orbit": {"type": "LEO", "altitude_km": 620, "INC": 90, "RAAN": 0, "TA": 0},
                          "pointing": {"type": "fixed", "dir": [0, 1, 0]}}}
}
```

Hay más ejemplos en `config/specs/` (SSO + GEO, formaciones, CW, Walker, TLE, órbitas aleatorias y modos de actitud en `ejemplo_actitud.json`).

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
│   └── test_geometry.jl            # Pruebas de las mallas de Geometry.jl (áreas, L_c, normales)
├── scripts/
│   ├── main.jl                     # Spec → JSON completo → simulación (propagación → actitud → fotometría → cámara → CSV)
│   ├── generate_config.jl          # Spec → JSON de escenario
│   ├── Kalman_RSO.jl               # Filtro de Kalman sobre las medidas LOS simuladas
│   ├── MonteCarlo.jl               # Campañas Monte Carlo sobre la configuración
│   ├── run_*_vscode.jl             # Lanzadores para VS Code (en el de simulación se elige la spec)
│   └── reduce_hipparcos_catalog.py # Reducción del catálogo Hipparcos
└── src/
    ├── Project.jl          # Módulo raíz: incluye todos los submódulos
    ├── Types.jl            # Tipos de datos (State, TrajectoryData, AttitudeData…)
    ├── Constants.jl        # Constantes físicas (μ, R, ω⊕, WGS84…)
    ├── Transformations.jl  # Rotaciones y conversiones de sistemas de referencia
    ├── Attitude/           # Actitud de los objetivos (ver abajo)
    │   ├── Attitude.jl       # Tipo base, registro de modos, lectura del bloque "attitude"
    │   ├── Quaternions.jl    # Cuaternión, Euler, DCM elementales y sorteo uniforme (T6)
    │   ├── Frames.jl         # LVLH, direcciones de apuntado y método TRIAD
    │   └── modes/            # spin, random, fixed, pointing
    ├── Environment.jl      # Lectura de configuración, SPICE, Sol, eclipses, LOS, exportación CSV
    ├── Propagator.jl       # Propagación orbital
    ├── Geometry.jl         # Geometría de los objetos (mallas, normales, áreas, coeficientes ópticos)
    ├── Photometry.jl       # Magnitud aparente (modelo de Phong por facetas)
    ├── Camera_sim.jl       # Generación de imagen: estrellas, PSF, fondo, ruido
    ├── Visualization.jl    # Visualización 3D de órbitas (GLMakie)
    └── Scenario/           # Generador de escenarios (ver abajo)
        ├── Scenario.jl       # Registros de plugins y autocarga
        ├── Builder.jl        # Valores por defecto, validación y escritura del JSON
        ├── Orbital_utils.jl  # Utilidades orbitales (SSO, LTAN, Kepler, r-v ↔ elementos)
        ├── Relative_utils.jl # Movimiento relativo (LVLH → elementos)
        ├── Orbits/           # keplerian, LEO, SSO, GEO, random_LEO
        ├── Patterns/         # single, formation, walker, tle
        ├── Layouts/          # along_track, element_offsets, cw_state, cw_ellipse, pco
        └── Pointing/         # fixed
```

### Generador de escenarios (`Scenario`)

Convierte una *spec* corta en el JSON completo que lee `Environment.load_config`. Los valores que la spec no indica se toman de `DEFAULT_TEMPLATE` y `DEFAULT_MEMBER_FIELDS` (`Builder.jl`), los mismos que lista `config/specs/plantilla.json`. Se basa en cuatro tipos de *plugin*, cada uno en su carpeta:

| Tipo | Qué define | Tipos disponibles |
|---|---|---|
| `Orbits/` | Elementos keplerianos de un objeto | `keplerian`, `LEO`, `SSO`, `GEO`, `random_LEO` |
| `Patterns/` | Cómo se generan los miembros de un grupo | `single`, `formation`, `walker`, `tle` |
| `Layouts/` | Disposición de una formación respecto al *chief* | `along_track`, `element_offsets`, `cw_state`, `cw_ellipse`, `pco` |
| `Pointing/` | Dirección de apuntado del observador | `fixed` |

Para añadir un tipo nuevo basta con crear un archivo en la carpeta correspondiente y registrarlo con `register_orbit!`, `register_pattern!`, `register_layout!` o `register_pointing!`. No hay que modificar nada más.

Cada grupo de objetivos debe definir su actitud con un bloque `attitude`, que se aplica a todos los miembros del grupo. Las claves antiguas `w_deg_s`/`w_dir` dan un error que muestra el bloque equivalente. En el modo `random`, cada miembro recibe su propia semilla a partir de la `seed` de la spec.

### Actitud de los objetivos (`Attitude`)

Cada objetivo del JSON del simulador lleva un bloque `attitude` obligatorio con un modo (`type`). El módulo lo valida y calcula la DCM cuerpo → ECI en cada instante:

| Modo | Parámetros | Comportamiento |
|---|---|---|
| `spin` | `w_deg_s`, `w_dir` (ejes cuerpo), orientación inicial opcional | Giro con ω constante alrededor de un eje fijo en el cuerpo |
| `random` | `w_range_deg_s` = [mín, máx], `seed` | Orientación inicial uniforme, eje uniforme y \|ω\| uniforme en el rango; luego como `spin` |
| `fixed` | `frame` (`ECI` o `LVLH`), orientación | Orientación constante respecto a ECI o a la órbita |
| `pointing` | `primary` y `secondary`: `body_axis` + `target` | Eje primario hacia `nadir`, `zenith`, `velocity`, `orbit_normal`, `sun` (con `-` para invertir) o un vector ECI; el secundario fija el giro (TRIAD) |

La orientación se da como cuaternión `q` o como `euler_deg` + `sequence`. Para añadir un modo basta con crear un archivo en `src/Attitude/modes/` y registrarlo con `register_attitude!`. Detalle de cada modo: `?parse_spin`, `?parse_random`, `?parse_fixed`, `?parse_pointing`.

## Modelos físicos

- **Propagación orbital:** gravedad central más armónicos zonales y teserales EGM2008 hasta grado y orden 4 (Montenbruck & Gill, 2000), integrada con `Vern9` (tolerancias 1e-11). Unidades: km, km/s; sistema inercial ECI (J2000).
- **Actitud de los objetivos:** convenciones del tema T6 (*Attitude Kinematics*). Cuaternión escalar primero, `C(q) = (q_η² − q_vᵀq_v) I + 2 q_v q_vᵀ − 2 q_η [q_v×]` como DCM referencia → cuerpo; Euler con secuencia `abc` y `C = C_c(θ3) C_b(θ2) C_a(θ1)`. El simulador usa la DCM cuerpo → ECI (`AttitudeData.DCM_b2i`, la transpuesta).
  - `spin`/`random`: solución exacta de la cinemática con ω constante, `C_ib(t) = C_ib(0) R(ê, ω t)`. No es dinámicamente consistente salvo que ê sea un eje principal de inercia (aún no hay dinámica con tensor de inercia).
  - `random`: orientación uniforme (cuaternión a partir de cuatro normales), eje uniforme en la esfera y |ω| uniforme en el rango.
  - LVLH (igual que en `Relative_utils.jl`): x radial, z normal al plano orbital (r × v), y = z × x.
  - `pointing`: método TRIAD con r, v (ECI) y la dirección del Sol (geocéntrica, SPICE). Error si las dos direcciones son paralelas en algún instante.
- **Entorno:** dirección del Sol con SPICE (DE440s), eclipses por la Tierra, línea de visión y cono de deslumbramiento solar.
- **Geometría de los objetos:** cuerpos triangulados (prismas, esferas, placas, cilindros, conos, troncos y fragmentos) con paneles opcionales y coeficientes ópticos aleatorios. Las dimensiones aleatorias siguen normales truncadas a valores positivos. Conos y troncos: solo superficie lateral exterior, sin bases ([PENDIENTE] confirmar si deben cerrarse). Fragmentos: envolvente convexa de puntos aleatorios escalada a la longitud característica del NASA Standard Breakup Model, L_c = (X + Y + Z)/3 (dimensiones en los ejes principales). Los paneles no tienen geometría propia: solo modifican los coeficientes ópticos medios.
- **Fotometría:** modelo de Phong por facetas (componentes difusa y especular) para obtener la magnitud aparente V.
- **Cámara:** catálogo Hipparcos dentro del campo de visión, PSF gaussiana, fondo, ruido Poisson, ruido de lectura y *bias*. Las medidas LOS se exportan con error. Por cada observador se escriben siempre `los_measurements.csv`, `state_measurements.csv` y `debris_positions.csv` (RSOs en el campo de visión con su posición en píxeles y su magnitud V en cada instante); los PNG solo con `image_generation_flag = true`.
- **Escenarios:**
  - Órbitas SSO con inclinación calculada por J2 y RAAN a partir de la LTAN.
  - Formaciones con condiciones iniciales de Clohessy-Wiltshire, válidas para un *chief* casi circular (e < 0.01) y separaciones ≪ SMA. Después se propagan con la dinámica completa.
  - Constelaciones Walker T/P/F (*delta* y *star*).
  - TLE propagados con SGP4 y convertidos de TEME a J2000 (IAU-76/FK5, sin EOP).

## Validación y resultados

- `test/test_attitude.jl` (41 comprobaciones superadas con Julia 1.13.2): matrices de Euler 3-2-1 y 3-1-3 y cuaternión frente a las expresiones del T6; el modo `spin` reproduce el antiguo `propagate_attitude` (diferencia < 1e-10) y cumple `dC_ib/dt = C_ib [ω×]`; media de las DCM aleatorias ≈ 0 (orientación uniforme); `fixed` en LVLH y `pointing` alinean los ejes pedidos. Se ejecuta con `julia --project=. test/test_attitude.jl`.
- `test/test_geometry.jl` (190 comprobaciones superadas con Julia 1.13.2): área de la esfera ≈ 4πr² (−1,2 % por ser un poliedro inscrito), área lateral del cono y del tronco frente a la analítica (< 1 %), L_c del fragmento igual a la pedida, ausencia de normales NaN y de dimensiones negativas en 24 000 sorteos, y generación de las 7 formas en todas las categorías. Se ejecuta con `julia --project=. test/test_geometry.jl`.
- Resto de la herramienta: [PENDIENTE]

## Créditos

- Autor de esta iteración: Brendan de la Oliva Yanguas.
- Autores de iteraciones anteriores: Ignacio Ricardo Ferrando Magdalena y Alejandro Ortiz Tello.
- Tutor: Ángel Luis Porras Hermoso.
- Máster Universitario en Sistemas Espaciales (MUSE), ETSIAE – Universidad Politécnica de Madrid.
