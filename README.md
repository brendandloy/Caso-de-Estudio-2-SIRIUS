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
- **Dinámica de actitud de los objetos.** No iniciado.
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

**1. Generar la configuración del escenario** a partir de una *spec* corta:

```bash
julia --project=. scripts/generate_config.jl config/specs/ejemplo_leo.json config/generated/ejemplo_leo.json
```

**2. Lanzar la simulación** con la configuración generada:

```bash
julia --project=. scripts/main.jl config/generated/ejemplo_leo.json output/
```

Sin argumentos, `main.jl` usa `config/config.json` y escribe en `output/`. Desde VS Code se puede ejecutar `scripts/run_simulation_vscode.jl`.

Ejemplo de *spec* (observador en LEO polar, objetivo en órbita kepleriana):

```json
{
  "scenario":   {"init_epoch": "2000-01-01T11:59:28", "Duration_hours": 0.2, "Data_time_step_seconds": 2},
  "objectives": {"deb":  {"orbit": {"type": "keplerian", "SMA": 8800, "ECC": 0.2, "INC": 90, "RAAN": 10, "AOP": 0, "TA": 357.5}, "w_deg_s": 0}},
  "observers":  {"cam1": {"orbit": {"type": "LEO", "altitude_km": 620, "INC": 90, "RAAN": 0, "TA": 0},
                          "pointing": {"type": "fixed", "dir": [0, 1, 0]}}}
}
```

Hay más ejemplos en `config/specs/` (SSO + GEO, formaciones, CW, Walker, TLE, órbitas aleatorias).

## Estructura del repositorio

```
├── config/
│   ├── specs/         # Specs cortas de escenario (entrada del generador)
│   ├── generated/     # Configuraciones completas generadas a partir de las specs
│   ├── tle/           # Archivos TLE de ejemplo
│   └── *.json         # Configuraciones heredadas (escritas a mano), Kalman y Monte Carlo
├── data/              # Kernels SPICE y catálogo Hipparcos
├── output/            # Resultados de simulación (ignorado por git)
├── scripts/
│   ├── main.jl                     # Simulación completa: propagación → fotometría → cámara → CSV
│   ├── generate_config.jl          # Spec → JSON de escenario
│   ├── Kalman_RSO.jl               # Filtro de Kalman sobre las medidas LOS simuladas
│   ├── MonteCarlo.jl               # Campañas Monte Carlo sobre la configuración
│   ├── run_*_vscode.jl             # Lanzadores para VS Code
│   └── reduce_hipparcos_catalog.py # Reducción del catálogo Hipparcos
└── src/
    ├── Project.jl          # Módulo raíz: incluye todos los submódulos
    ├── Types.jl            # Tipos de datos (State, TrajectoryData, AttitudeData…)
    ├── Constants.jl        # Constantes físicas (μ, R, ω⊕, WGS84…)
    ├── Transformations.jl  # Rotaciones y conversiones de sistemas de referencia
    ├── Environment.jl      # Lectura de configuración, SPICE, Sol, eclipses, LOS, exportación CSV
    ├── Propagator.jl       # Propagación orbital y de actitud
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

Convierte una *spec* corta en el JSON completo que lee `Environment.load_config`. Se basa en cuatro tipos de *plugin*, cada uno en su carpeta:

| Tipo | Qué define | Tipos disponibles |
|---|---|---|
| `Orbits/` | Elementos keplerianos de un objeto | `keplerian`, `LEO`, `SSO`, `GEO`, `random_LEO` |
| `Patterns/` | Cómo se generan los miembros de un grupo | `single`, `formation`, `walker`, `tle` |
| `Layouts/` | Disposición de una formación respecto al *chief* | `along_track`, `element_offsets`, `cw_state`, `cw_ellipse`, `pco` |
| `Pointing/` | Dirección de apuntado del observador | `fixed` |

Para añadir un tipo nuevo basta con crear un archivo en la carpeta correspondiente y registrarlo con `register_orbit!`, `register_pattern!`, `register_layout!` o `register_pointing!`. No hay que modificar nada más.

## Modelos físicos

- **Propagación orbital:** gravedad central más armónicos zonales y teserales EGM2008 hasta grado y orden 4 (Montenbruck & Gill, 2000), integrada con `Vern9` (tolerancias 1e-11). Unidades: km, km/s; sistema inercial ECI (J2000).
- **Actitud:** rotación a velocidad angular constante (módulo `w_deg_s`, eje `w_dir`) sobre una DCM cuerpo → inercial.
- **Entorno:** dirección del Sol con SPICE (DE440s), eclipses por la Tierra, línea de visión y cono de deslumbramiento solar.
- **Geometría de los objetos:** cuerpos triangulados (prismas, esferas, placas, cilindros, conos, troncos y fragmentos) con paneles opcionales y coeficientes ópticos aleatorios.
- **Fotometría:** modelo de Phong por facetas (componentes difusa y especular) para obtener la magnitud aparente V.
- **Cámara:** catálogo Hipparcos dentro del campo de visión, PSF gaussiana, fondo, ruido Poisson, ruido de lectura y *bias*. Las medidas LOS se exportan con error.
- **Escenarios:**
  - Órbitas SSO con inclinación calculada por J2 y RAAN a partir de la LTAN.
  - Formaciones con condiciones iniciales de Clohessy-Wiltshire, válidas para un *chief* casi circular (e < 0.01) y separaciones ≪ SMA. Después se propagan con la dinámica completa.
  - Constelaciones Walker T/P/F (*delta* y *star*).
  - TLE propagados con SGP4 y convertidos de TEME a J2000 (IAU-76/FK5, sin EOP).

## Validación y resultados

[PENDIENTE]

## Créditos

- Autor de esta iteración: Brendan de la Oliva Yanguas.
- Autores de iteraciones anteriores: Ignacio Ricardo Ferrando Magdalena y Alejandro Ortiz Tello.
- Tutor: Ángel Luis Porras Hermoso.
- Máster Universitario en Sistemas Espaciales (MUSE), ETSIAE – Universidad Politécnica de Madrid.
