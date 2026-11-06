# CLAUDE.md — SIRIUS (Caso de Estudio 2)

Simulador de imágenes de objetos residentes en el espacio (RSOs) captados por cámaras, en Julia. Esta iteración añade dinámica de actitud, mejores escenarios y más realismo en la magnitud de los objetos. Contexto completo, modelos físicos y estructura: `README.md`.

- Julia 1.13 (versiones fijadas en `Manifest.toml`). Paquete `SIRIUS` (`src/SIRIUS.jl`).
- Unidades del código: km, km/s, s, rad (grados solo en el JSON y donde el nombre lo indica, p. ej. `_deg`). Sistema de referencia: ECI (J2000) salvo que se indique.
- Estado actual: Kalman (`scripts/Kalman_RSO.jl`) y Monte Carlo (`scripts/MonteCarlo.jl`) no se usan en esta iteración; no los modifiques salvo que lo pida.

## Comandos

- Todas las pruebas: `julia --project=. test/runtests.jl` (o `] test`). Una: `julia --project=. test/test_<módulo>.jl`.
- Simulación: `julia --project=. scripts/main.jl config/specs/<spec>.json`.
- Rendimiento: `julia --project=. scripts/benchmark.jl [<spec>]` (tiempo y memoria, sin compilación).
- Dependencias: `] instantiate`; tras renombrar el paquete o cambiar dependencias, `] resolve`.

## Cómo escribir código

### Organización
- `src/`: toda la lógica, sin duplicar. `scripts/`: solo cargan, llaman a funciones de `src/` y guardan resultados; no definen funciones salvo helpers mínimos.
- Antes de escribir una función, busca si ya existe (p. ej. lectores del JSON en `JSON_utils.jl`, conversiones en `Transformations.jl`, constantes en `Constants.jl`, parámetros del programa en `Parameters.jl`).
- Constantes físicas solo en `Constants.jl`; valores que fijamos los desarrolladores, en `Parameters.jl`; lo que elige el usuario, en el JSON.
- Rutas siempre relativas a la raíz (`DATA_DIR`, `@__DIR__`), nunca absolutas.
- Las claves del JSON son la interfaz con el usuario: no se renombran sin avisar, porque rompen las specs.
- Ampliaciones con el patrón de registro existente (`register_attitude!`, `register_shape!`, `register_orbit!`...): un archivo nuevo por modo, forma u órbita.

### Funciones
- Una sola tarea por función (si al describirla hace falta "y además", se divide) y un solo nivel de abstracción: las funciones de alto nivel llaman a pasos con nombre, sin bucles de bajo nivel mezclados.
- Orientativo: más de ~30 líneas de código es señal de que conviene dividir, salvo que la geometría o el algoritmo lo pidan (mallas, recursiones).
- Pocos argumentos posicionales; si son muchos, agrúpalos (struct o NamedTuple) o usa keywords.
- `struct` cuando los campos son fijos y se conocen al escribir el código (`DebrisModel`), con campos de tipo concreto; `Dict` cuando las claves solo se conocen al ejecutar (trayectorias por nombre de objeto, registros de plugins, JSON recién leído). Nunca un `Dict{Symbol,Any}` con claves fijas.
- Sin código muerto ni comentado: git conserva el historial.

### Nombres
- `snake_case` para funciones y variables; `CamelCase` para módulos y tipos; `UPPER_SNAKE_CASE` para constantes (`MU_EARTH`, `R_EARTH`).
- Booleanos con `is_`/`has_` (`is_visible`, `has_los`); `needs_` también vale.
- Las siglas siguen el estilo del nombre: minúscula en snake_case (`has_los`, `is_pole_in_fov`), mayúscula en CamelCase (`LOSMeasurement`, `JSONUtils`) y en constantes (`JD_J2000`).
- Nombres que expliquen la intención; nada de abreviaturas crípticas (`sf`, `b`). Los símbolos de las ecuaciones (μ, ω, ρ, θ) sí, si coinciden con la documentación.
- Unidades en el nombre cuando no sean las del código (`w_deg_s`, `dims_m`, `area_km2`).

### Documentación
- Cabecera al inicio de cada módulo y script: qué hace, funciones o archivos que contiene, dependencias y cualquier convenio relevante.
- `function ... end`: docstring justo encima de `function`, sin líneas en blanco (para `?nombre`), con:
  1. Firma y descripción en una o dos frases.
  2. `Por qué existe`: para qué se usa y por qué es una función aparte (una o dos líneas; se omite si es evidente).
  3. `# Argumentos` (tipo, unidades, sistema de referencia), `# Devuelve`, `# Excepciones`.
  4. `# Ejemplo` solo en funciones públicas o con convenios que se puedan confundir (cuaterniones, DCM, signos); máximo 2–3 líneas.
  5. Si aplica: ecuación del modelo y su referencia (p. ej. M&G (2000) ec. 3.29; T6, diap. 49).
- Sintaxis compacta (`f(x) = ...`): comentario de una o dos líneas si es un helper interno y obvio; docstring corto (firma, una o dos frases, unidades) si se usa desde otro archivo o lleva un convenio.
- Comentarios de una o dos líneas que expliquen el porqué, las unidades y los convenios; no repitan lo que el código ya dice.

### Pruebas y reproducibilidad
- Cada módulo de `src/` tiene su `test/test_<módulo>.jl` con casos analíticos calculados a mano (nunca datos reales); se añade a `TEST_FILES` de `runtests.jl`.
- Una corrección de bug va con el test que lo detecta (primero el test, que debe fallar; después el arreglo).
- Todo lo aleatorio con semilla (`get_or_draw_seed`): misma spec y semilla, mismo resultado.
- Optimizar solo tras medir (`scripts/benchmark.jl`, `Profile`) y sin perder claridad. `@inbounds` solo con longitudes comprobadas.

## Cómo trabajar

- Cambios mínimos y localizados; nada de reescrituras completas sin pedirlo. Lee siempre la versión actual de un archivo antes de editarlo.
- Conserva los finales de línea de cada archivo (la mayoría usan CRLF de Windows).
- Trabajo en ramas (`nombre_tarea`), nunca directamente en `main`, y por tandas: una tanda = un cambio que se resume en una frase = un commit.
- Cada tanda se verifica antes del commit: `runtests.jl` completo y, si el cambio no debe alterar resultados, comparación de `output/ejemplo_actitud/` antes y después (`rso_states.csv`, `debris_positions.csv` y columnas verdaderas de `state_measurements.csv` y `los_measurements.csv`; las columnas `*_meas` llevan ruido).
- No incluir en los commits `output/` ni los `config/generated/` de specs sin semilla (cambian en cada ejecución).
## README y CHANGELOG

- Se actualizan en la misma tanda en que se crea, renombra o elimina un módulo, script, plugin o función pública; cambia un modelo físico o un supuesto; cambian dependencias, instalación o forma de ejecutar; o se obtiene un resultado o validación. Si me dices que has hecho commits por tu cuenta, revisa el historial posterior a la última entrada del CHANGELOG y documenta lo que falte.
- Cambios mínimos en las secciones afectadas; no copiar docstrings (resumir y remitir a `?nombre_funcion`); unidades y sistemas de referencia explícitos; no inventar resultados ni funcionalidades; lo que falte, `[PENDIENTE]`.
- README (español): 1) título, descripción en 2–3 líneas y estado; 2) contexto: casos de estudio anteriores y objetivos de esta iteración con su estado (En curso / Completado / No iniciado); 3) requisitos e instalación; 4) uso rápido con ejemplo ejecutable; 5) estructura del repositorio (árbol comentado, una línea por archivo o módulo); 6) modelos físicos: supuestos, sistemas de referencia y unidades; 7) validación y resultados; 8) créditos: autor (Brendan de la Oliva Yanguas), autores anteriores (Ignacio Ricardo Ferrando Magdalena, Alejandro Ortiz Tello), tutor (Ángel Luis Porras Hermoso), MUSE – ETSIAE – UPM. Sin imágenes, licencia, limitaciones ni rama de trabajo.
- CHANGELOG: formato Keep a Changelog, una sección por fecha (lo más reciente arriba) con Añadido / Cambiado / Corregido / Eliminado. Solo lo que note quien usa o continúa la herramienta, no cada commit.

## Commits

Título y descripción para GitHub Desktop:

```
<Verbo en 3.ª persona del presente> <qué cambia>        (≤ 72 caracteres, sin punto)

<Párrafo de 2–3 líneas: qué cambia y por qué; si los resultados cambian o no.>

Añadido:
- <archivo/función>: <qué hace>.

Cambiado:
- ...

Corregido:
- ...

Eliminado:
- ...

Verificación:
- <runtests.jl: número de comprobaciones>.
- <comparación de resultados, simulaciones o specs probadas>.

README y CHANGELOG actualizados.
```

- Solo las secciones que apliquen, en ese orden; líneas de ~72 caracteres. `Verificación:` va siempre: qué se ha comprobado y con qué resultado.
- Ejemplos de título: "Corrige el instante de más en propagate cuando T no es múltiplo de Δt", "Añade bloque shape y modelo de población para los objetivos".
