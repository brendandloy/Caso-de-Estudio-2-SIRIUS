# Changelog

Cambios relevantes de SIRIUS en el Caso de Estudio 2. Formato basado en [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/); lo más reciente, arriba.

## [Sin publicar]

### 2026-11-06
#### Añadido
- `scripts/benchmark.jl`: mide con `@time` la propagación y la simulación completa de una spec (sin la compilación), para comparar el rendimiento antes y después de un cambio.
- `Geometry.DebrisModel`: tipo con los campos del modelo fotométrico de un objetivo (ver `?DebrisModel`).
- `src/Visibility.jl` (módulo `Visibility`): eclipse por la Tierra, línea de visión y cono de deslumbramiento solar, con el modelo y las ecuaciones en sus docstrings (`?eclipseEarthECI`, `?losECI`, `?sunGlareConeECI`).
- `test/test_visibility.jl`: pruebas geométricas calculadas a mano de las tres funciones.
- `src/JSON_utils.jl` (módulo `JSONUtils`): lectores y validadores comunes de los bloques del JSON (`todict`, `check_keys`, `getnum`, `getvec`, `getunit`), que antes estaban en Attitude y los importaba Geometry, y `get_or_draw_seed` (ver `?get_or_draw_seed`), que lee la semilla de un bloque o sortea una nueva y la imprime. Se carga antes que Attitude; `test_attitude.jl`, `test_geometry.jl` y `generate_config.jl` lo incluyen.
- `test/test_propagator.jl`: comprueba que la trayectoria propagada tiene los mismos instantes que la malla 0:Δt:T de `scripts/main.jl`.
- `Constants.DATA_DIR`: ruta única a la carpeta `data/` de la raíz; la usan `scripts/main.jl` (kernels SPICE) y `Camera_sim.jl` (catálogo Hipparcos) en lugar de escribirla cada uno.

#### Cambiado
- `generate_debris` devuelve un `DebrisModel` en lugar de un `Dict{Symbol,Any}`: los campos se leen como `debris.normals` en vez de `debris[:normals]`. `compute_apparent_magnitude` y `_phong_model` lo reciben entero en lugar de cinco valores sueltos.
- `_phong_model` suma el flujo faceta a faceta en un bucle, sin arrays temporales; misma ecuación y mismo resultado (la dirección especular ya era unitaria, así que se quita su normalización).
- `accel_earth_perturbations`: el grado máximo por defecto es el de la tabla de coeficientes (`size(C, 1) - 1`), en lugar de 4 fijo.
- `eclipseEarthECI`, `losECI` y `sunGlareConeECI` pasan de `Environment.jl` a `Visibility.jl`; Photometry las toma de ahí. El resultado no cambia, solo lo que devuelven: `eclipseEarthECI` → `(en_sombra, region)`, `losECI` → `Bool`, `sunGlareConeECI` → `(deslumbrado, angulo_deg)`. `Environment.jl` queda con la lectura de la configuración, SPICE, la dirección del Sol, el error de navegación y la exportación de CSV.
- Un único `check_keys(b, allowed; who, computed=())` para Attitude, Geometry y Scenario (antes había uno en Attitude con `who` posicional y otro en `Orbital_utils.jl`). Lista todas las claves no admitidas a la vez; el mensaje de las órbitas cambia ligeramente de formato.
- La semilla de la actitud `random`, de la forma de los objetivos y de la spec se lee con `get_or_draw_seed`; los avisos de semilla sorteada tienen ahora el mismo formato en los tres casos.
- Las configuraciones heredadas escritas a mano (`config/config*.json`) pasan a `config/legacy/`. `mc_config.json` apunta a `legacy/config _cross_e.json`.
- `axis_angle_rotation` (fórmula de Rodrigues) pasa de `Transformations.jl` a `test/test_attitude.jl`: solo la usan esas pruebas, como solución de referencia independiente de los cuaterniones.
- Docstring de `accel_earth_perturbations` colocado encima de la función (antes estaba dentro y `?accel_earth_perturbations` no lo mostraba) y con las unidades reales (km, km/s²).

#### Corregido
- `propagate` devolvía un instante de más (el final T) cuando la duración no era múltiplo del paso de salida, y la trayectoria quedaba con una muestra más que la dirección del Sol; la fotometría leía entonces fuera de `sun_dir` sin dar error. Ahora la salida es exactamente 0:Δt:T. Las simulaciones con T múltiplo de Δt no cambian.

#### Eliminado
- Constante `DCM_I` sin uso en `Photometry.jl`.
- Opción `return_details` de las tres funciones (tupla de diagnóstico con el motivo del bloqueo y las distancias intermedias). No se usaba en ningún sitio; se puede recuperar del historial de git anterior a este cambio.
- Dependencias sin uso `AstroTime` y `QHull` (la envolvente convexa es propia desde el 2026-11-01).
- Código sin uso: `wrap_pi` y `add_bias!` (`Camera_sim.jl`), el tipo `StationTrack` (`Types.jl`), `Constants.data_path` (apuntaba a `src/data`, que no existe) y `phi_sun` comentado en `Constants.jl` (la irradiancia se lee de `photometry_config.phi_sun`).

### 2026-11-04
#### Añadido
- Actitud de los observadores: bloque `attitude` obligatorio en cada observador (spec y JSON del simulador), con los mismos modos que los objetivos. Describe la actitud del satélite.
- Modo de actitud `tracking` (`src/Attitude/modes/tracking.jl`, solo observadores): el boresight de la cámara apunta en cada instante a un objetivo (`objective`, nombre del miembro en el JSON) y `roll` (`body_axis` + `direction`) fija el giro del satélite alrededor de ese eje.
- Bloque `camera_mount` de los observadores (`src/Attitude/Mount.jl`): `boresight` y `up` en ejes cuerpo definen el montaje de la cámara, C_cb. Sin bloque se usa boresight = +z cuerpo y arriba = −y cuerpo (C_cb = I) y se avisa; el generador lo escribe en el JSON.
- `bind_attitude`: paso opcional de los modos que necesitan datos externos (montaje de la cámara, trayectorias de los objetivos); `camera_dcm_history` calcula C_ci = C_cb · C_bi.
- Pruebas de `camera_mount`, `tracking` y orientación de la imagen en `test/test_attitude.jl`.

#### Cambiado
- La cámara recibe su orientación C_ci (ECI → cámara) en cada instante en lugar de una dirección fija. Proyección pinhole con ejes cámara dextrógiros (x derecha, y abajo, z boresight), con la misma focal que antes; con actitud variable el campo estelar se recalcula en cada imagen. Corrige el espejo de las imágenes anteriores: con el norte arriba, el este sale ahora a la izquierda, así que los píxeles de `debris_positions.csv` cambian. Con actitud variable, `observed_star_catalog.csv` añade la columna `image_name`.
- Modo `pointing`: la clave `target` pasa a llamarse `direction` (también en `Frames.jl`: `parse_direction`, `direction_vector`).
- Specs de ejemplo y `plantilla.json`: el apuntado fijo se escribe como `attitude` `pointing` hacia el mismo vector ECI, con el norte celeste arriba, y `camera_mount` explícito. `ejemplo_actitud.json` añade una cámara en `tracking`.

#### Eliminado
- Bloque `pointing` de los observadores, `obs_dir` en el JSON del simulador y los plugins de apuntado del generador (`src/Scenario/Pointing/`, `register_pointing!`).
- Mensaje específico para las claves antiguas `w_deg_s`/`w_dir` y `legacy_to_spin`: ahora dan el error general de bloque `attitude` ausente.

### 2026-11-01
#### Añadido
- Bloque `shape` obligatorio en cada objetivo (spec y JSON del simulador): `type` (`sphere`, `prism`, `plate`, `cylinder`, `cone`, `trunk`, `fragment` o `random`) y, opcionalmente, todas sus medidas en metros. Sin medidas se sortean; con `random` se sortean también la clase de objeto y la forma. El generador escribe el bloque ya resuelto, con una `seed` por miembro.
- Modelo de población para las formas sorteadas (`src/Geometry/Population.jl`), con datos del informe de entorno espacial de ESA de 2026: clases de objeto de la tabla 3.1, ley de potencias por tramos de MASTER 2026 y del NASA SBM para los fragmentos, e histogramas de área de satélites y etapas.
- Parámetros maestros de la forma en `src/Parameters.jl`, entre ellos `SHAPE_LC_MIN_M = 0.10` m. Con medidas concretas por debajo de ese valor se avisa de que el objetivo probablemente no se verá.
- Ejes cuerpo y longitud característica documentados para cada forma (cabecera de `src/Geometry/Geometry.jl` y README).
- Registro de formas (`register_shape!`): para añadir una basta con un archivo en `src/Geometry/shapes/`.

#### Cambiado
- `Geometry` pasa a la carpeta `src/Geometry/` (`Geometry.jl`, mallas, población y una forma por archivo), como `Attitude` y `Scenario`, y se carga antes de `Scenario`, que lo usa para resolver las formas. El antiguo `src/Geometry.jl` se elimina.
- `generate_debris(spec::ShapeSpec)` sustituye a `generate_debris(sat_size, sub_cat; deploy_panels, shape_type)`. Malla del fragmento, coeficientes ópticos y `alpha` usan la semilla del objetivo: la misma spec da siempre el mismo objeto.
- Tronco de cono: el radio superior es una medida (`top_radius_m`) en lugar de sortearse dentro de la malla.
- Prisma sorteado: el eje z es la dirección longitudinal (z ≥ x ≥ y), de modo que las caras ±z (tapas, con su propio material) son los extremos.
- Fragmento: la malla se expresa en sus ejes principales (x la dimensión mayor) y la envolvente convexa se calcula sin QHull (comparada con `scipy.spatial.ConvexHull`). QHull sigue en `Project.toml`, pero ya no se usa.
- Specs de ejemplo y `plantilla.json` con su bloque `shape`. `ejemplo_leo.json` conserva la caja de 4,1 × 4,1 × 12,8 m que se simulaba antes.
- `test/test_geometry.jl` reescrito para el nuevo módulo (602 comprobaciones).

#### Eliminado
- `sat_size`, `sub_cat` y las funciones `_dims_*` (catálogo de medidas heredado); `deploy_panels` y la caja fija de 4,1 × 4,1 × 12,8 m escrita a mano en `generate_debris`; `generateRandomShape`. Los paneles solares se modelarán más adelante.

### 2026-10-30
#### Añadido
- `src/Parameters.jl`: parámetros maestros del programa (configuración que no elige el usuario). Primer parámetro: `RANDOM_W_RANGE_DEG_S = (0.5, 5.0)`, rango de |ω| del modo de actitud `random`.
- `dcm_to_quat` (método de Shepperd), `omega_matrix` y `quat_step` en `Attitude/Quaternions.jl`.

#### Cambiado
- Modos `spin` y `random`: la actitud se propaga como cuaternión con la cinemática del T6 (`q̇ = ½ Ω(ω) q`), paso a paso con la exponencial exacta para ω constante y renormalización. Los resultados coinciden con la versión anterior (diferencia del orden del redondeo).
- Modo `random`: el rango de |ω| ya no se da en el JSON (`w_range_deg_s` da error); se toma de `RANDOM_W_RANGE_DEG_S`. Ejemplo `ejemplo_actitud.json` actualizado.

### 2026-10-29
#### Añadido
- Pruebas `test/test_geometry.jl`.

#### Cambiado
- Carpeta de resultados por defecto de `main.jl`/`run_spec`: `output/<nombre de la spec>/` en lugar de `output/`. En `run_simulation_vscode.jl`, `output = nothing` usa esa carpeta.

#### Corregido
- Fragmento irregular (forma 7): daba error siempre, porque llamaba a una función `qhull` que QHull.jl no exporta; ahora usa `chull`. Además, la longitud característica se usaba como si fuera un volumen (un fragmento de 0.1 m salía con ≈ 29 m de radio); ahora la malla se escala para que L_c = (X + Y + Z)/3 sea la pedida.
- Conos y troncos con el área duplicada: las caras interiores repetían las exteriores. Ahora son ≈ 0.75 mag más débiles, el valor correcto para su tamaño. Se mantienen sin bases (superficies abiertas).
- Dimensiones aleatorias negativas en prismas, conos y troncos (≈ 2 % de los conos pequeños): ahora se sortean de una normal truncada a valores positivos.
- Esfera: los polos ya no se representan con anillos de puntos casi coincidentes, sino con un vértice y un abanico de triángulos.
- `calculateNormals` descarta, con un aviso, los triángulos degenerados (área nula) en lugar de devolver normales NaN.

### 2026-10-27
#### Añadido
- `config/specs/plantilla.json`: spec con todas las opciones del escenario (`scenario`, `member_fields`, `seed`) y su valor por defecto, para copiarla al crear un escenario nuevo.
- El generador rechaza claves desconocidas en la spec, en `scenario` (incluidos `camera_config` y `photometry_config`), en el `template` y en `member_fields`.

#### Cambiado
- `debris_positions.csv` se escribe siempre, aunque no se generen imágenes; `image_generation_flag` solo controla los PNG.
- README: tabla con las tres formas de ejecutar la simulación y aviso de que el REPL abierto no recoge cambios en `src/`.

#### Eliminado
- Comprobación específica de `w_deg_s`/`w_dir` en `member_fields` (la cubre la comprobación general de claves).

### 2026-10-26
#### Añadido
- Módulo `Attitude` (`src/Attitude/`): actitud de los objetivos con cuatro modos, `spin` (ω constante), `random` (orientación, eje y |ω| aleatorios con semilla), `fixed` (orientación constante respecto a ECI o LVLH) y `pointing` (eje primario y secundario hacia nadir, Sol, velocidad, normal orbital o un vector ECI; método TRIAD). Orientaciones como cuaternión o ángulos de Euler, con las convenciones del tema T6.
- Bloque `attitude` obligatorio en los objetivos de las specs del generador. En `random`, cada miembro recibe su propia semilla a partir de la `seed` de la spec.
- Ejemplo `config/specs/ejemplo_actitud.json` con los cuatro modos.
- Pruebas `test/test_attitude.jl`.

#### Cambiado
- `scripts/main.jl` recibe una spec (obligatoria) y la carpeta de resultados (opcional, `output/` por defecto): genera el JSON completo en `config/generated/` y simula con él. Ya no tiene valor por defecto `config/config.json`. `run_simulation(json, output)` sigue disponible para JSON completos (lo usa `MonteCarlo.jl`).
- `run_simulation_vscode.jl`: la spec a simular se elige con la variable `spec`.
- El módulo `Scenario` forma parte del módulo `Project`.
- El simulador y el generador exigen el bloque `attitude` en cada objetivo; `w_deg_s`/`w_dir` dan un error con el bloque equivalente. Los JSON heredados de `config/` siguen en el formato antiguo.
- Specs de ejemplo migradas al bloque `attitude` y `config/generated/` regenerado.
- La spec ya no admite `w_deg_s`/`w_dir` en `member_fields`: la actitud se define en cada grupo de objetivos.

#### Corregido
- El eje de giro `w_dir` se normaliza; antes un eje no unitario (p. ej. `[1, 1, 0]`) daba una rotación incorrecta.
- `load_config` lee el JSON en memoria en vez de mapearlo (mmap). En Windows, el mapeo impedía sobrescribir el JSON de `config/generated/` mientras siguiera abierta la sesión de Julia que lo había leído (error `Invalid argument`).

#### Eliminado
- `Propagator.propagate_attitude`, sustituido por `Attitude.compute_attitude`.

### 2026-10-23
#### Añadido
- Tipo de órbita `random_LEO`: LEO circular aleatoria. Admite rangos o valores fijos para la altitud, INC, RAAN y TA, con INC isótropa (cos INC uniforme).
- Semilla `seed` en la spec para que los escenarios aleatorios sean reproducibles.
- Ejemplo `config/specs/ejemplo_random.json`.

### 2026-10-22
#### Añadido
- Patrón `tle`: un miembro por cada TLE de un archivo, propagado con SGP4 hasta `init_epoch` y convertido de TEME a J2000. Admite selección por NORAD ID, por nombre o por rango.
- Dependencias `SatelliteToolboxTle`, `SatelliteToolboxSgp4` y `SatelliteToolboxTransformations`.
- Patrón `formation` con las disposiciones `along_track`, `element_offsets`, `cw_state`, `cw_ellipse` y `pco` (Clohessy-Wiltshire).
- Patrón `walker`: constelaciones Walker T/P/F (*delta* y *star*).
- `Relative_utils.jl`: paso de estado relativo LVLH a elementos keplerianos, con avisos de validez de CW.
- `CLAUDE.md` con las normas de documentación del código.
- Ejemplos `ejemplo_tle`, `ejemplo_formacion`, `ejemplo_cw` y `ejemplo_walker`.

### 2026-10-18
#### Añadido
- Tipos de órbita `SSO` (INC por J2 y RAAN a partir de la LTAN) y `GEO` (a partir de la longitud subsatélite).
- `Orbital_utils.jl`: utilidades orbitales compartidas (SSO, LTAN, ERA, Kepler, conversión r-v ↔ elementos).
- Ejemplo `ejemplo_sso_geo`.

#### Cambiado
- Refactorización de `Builder.jl` y mejora de los comentarios y docstrings de todo el módulo `Scenario`.

### 2026-10-17
#### Añadido
- Módulo `Scenario`: generador que convierte una spec corta en el JSON completo de escenario, con arquitectura de *plugins* (órbitas, patrones, disposiciones y apuntado).
- Órbitas `keplerian` y `LEO`, patrón `single` y apuntado `fixed`.
- Script `scripts/generate_config.jl` y ejemplo `ejemplo_leo`.

## [0.1.0] – 2026-10-13 – Versión heredada

Código de las iteraciones anteriores de SIRIUS (CE2 de I. R. Ferrando Magdalena y A. Ortiz Tello; CE3 de I. R. Ferrando Magdalena), subido como punto de partida:

- Propagación orbital (EGM2008 4×4) y de actitud (rotación uniforme).
- Entorno: SPICE, eclipses, línea de visión y deslumbramiento solar.
- Geometría de los objetos y fotometría de Phong.
- Simulación de la cámara: catálogo Hipparcos, PSF, fondo y ruido.
- Scripts de filtro de Kalman y Monte Carlo, y configuraciones escritas a mano.
