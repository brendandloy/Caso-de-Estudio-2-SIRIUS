# Changelog

Cambios relevantes de SIRIUS en el Caso de Estudio 2. Formato basado en [Keep a Changelog](https://keepachangelog.com/es-ES/1.1.0/); lo más reciente, arriba.

## [Sin publicar]

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
