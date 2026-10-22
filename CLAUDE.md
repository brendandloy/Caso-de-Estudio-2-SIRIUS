# CLAUDE.md - Proyecto SIRIUS (Caso de Estudio 2)

## Información del Proyecto
- **Título:** Desarrollo de herramienta para la simulación de imágenes de objetos residentes en el espacio D3 (SIRUS).
- **Descripción:** Continuación de la herramienta SIRIUS, diseñada para simular imágenes de RSOs captadas a través de cámaras. Esta iteración añade módulos de dinámica de actitud, mejoras en escenarios e incremento del realismo en la magnitud relativa de los objetos.

## Stack Tecnológico
- **Lenguaje:** Julia

---

## Directrices de Desarrollo y Estilo de Código

### 1. Documentación de Funciones (Docstrings)
Al generar funciones en Julia, sigue estrictamente estas reglas de documentación:

- **Sintaxis Estándar (`function ... end`):** Coloca un docstring formateado **justo encima** de la línea `function` (sin líneas en blanco intermedias) para que el REPL de Julia lo registre adecuadamente (accesible con `?nombre_funcion`).
  - **Contenido obligatorio:** Descripción de la función, argumentos de entrada, valor de retorno, excepciones/errores posibles y un ejemplo de uso si ayuda a la comprensión.
- **Sintaxis Compacta (`f(x) = x + 1`):** No uses docstrings. Coloca únicamente un comentario explicativo normal encima de la función.

### 2. Documentación de Módulos y Scripts
- **Módulos:** Al inicio del módulo, incluye una descripción detallada de lo que hace, los scripts que contiene y cualquier información relevante.
- **Scripts:** Al inicio de cada script, incluye una descripción de su propósito, las funciones que contiene y cualquier información relevante.

### 3. Comentarios en el Código
- Genera comentarios explicativos breves a lo largo de todo el script (en funciones, variables, bloques lógicos, etc.).
- Los comentarios deben tener una longitud máxima de **una o dos líneas**.

---

## Comandos Frecuentes de Julia (Referencia rápida)
- Iniciar REPL: `julia`
- Activar entorno del proyecto: `julia --project=.` o `] activate .` dentro del REPL.
- Ejecutar un script: `julia --project=. ruta/al/script.jl`
- Ejecutar tests (si aplican): `] test`
