# Changelog

Cambios relevantes de SIRIUS en el Caso de Estudio 2. Formato basado en [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/); lo más reciente, arriba.

## [Sin publicar]

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
