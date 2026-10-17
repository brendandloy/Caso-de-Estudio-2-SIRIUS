# Se incluye desde Scenario.jl (comparte su namespace).

const ELEMENT_KEYS = ("SMA", "ECC", "INC", "RAAN", "AOP", "TA")

# Valores por defecto (los de tus config_*.json). Se pisan con "template" y/o "scenario".
const DEFAULT_TEMPLATE = Dict{String,Any}(
    "Duration_hours"                   => 0.2,
    "Data_time_step_seconds"           => 2,
    "init_epoch"                       => "2000-01-01T11:59:28",
    "navigation_error_flag"            => true,
    "graphic_display_flag"             => false,
    "image_generation_flag"            => false,
    "ignore_fov_for_measurements_flag" => false,
    "star_csv_flag"                    => false,
    "camera_config" => Dict{String,Any}(
        "Vmag_lim"          => 10.0,
        "fov_deg"           => 27,
        "resolution_x"      => 1024,
        "resolution_y"      => 1024,
        "mag0"              => 8.0,
        "flux0"             => 20000.0,
        "fwhm_deg"          => 0.0659,
        "background_level"  => 500,
        "sigma_read_noise"  => 100.0,
        "upper_percentile"  => 99.9,
        "lower_percentage"  => 0.0,
        "los_sigma_arcsec"  => 42.43
    ),
    "photometry_config" => Dict{String,Any}(
        "phi_sun"         => 455.0,
        "glare_angle_deg" => 15.0,
        "Vref_phong"      => -10.0,
        "Fref_phong"      => 0.000001
    ),
    "stations_llh" => Dict{String,Any}()
)

# Valores por defecto de los campos de las camaras y objetivos.
const DEFAULT_MEMBER_FIELDS = Dict{String,Any}(
    "position_sigma_km"   => [0.015, 0.015, 0.015],
    "velocity_sigma_km_s" => [0.00002, 0.00002, 0.00002],
    "w_deg_s"             => 0.0,
    "w_dir"               => [1, 0, 0]
)

# -----------------------------------------------------------------------------
# JSON: lectura, fusión y escritura (con orden de claves estable)
# -----------------------------------------------------------------------------

function read_json(path::AbstractString)
    """
    Lee un archivo en formato JSON desde una ruta dada y lo convierte a un diccionario de Julia.

    # Argumentos
    - `path::AbstractString`: Ruta relativa o absoluta del archivo JSON a leer.

    # Devuelve
    - `Dict{String,Any}`: Diccionario con la estructura de datos leída del archivo JSON.

    # Excepciones
    - `ErrorException`: Se lanza si el archivo no existe en la ruta especificada.
    - `ErrorException`: Se lanza si ocurre un error al intentar interpretar el contenido del JSON.
    """

    # Verifica que el archivo exista en el disco; si no, lanza un error
    isfile(path) || error("File not found: $path")
    try
        # Lee el contenido del archivo como cadena de texto y lo parsea a un Dict{String,Any}
        return JSON3.read(read(path, String), Dict{String,Any})
    catch err
        # Captura cualquier fallo de sintaxis en el JSON y lanza un error explicativo
        error("Could not parse '$path' (revisa llaves sobrantes al final): " * sprint(showerror, err))
    end
end

function deep_merge!(base::AbstractDict, over::AbstractDict)
    """
    Combina recursivamente dos diccionarios, sobrescribiendo o fusionando sus contenidos en profundidad.

    # Argumentos
    - `base::AbstractDict`: Diccionario base que será modificado con los nuevos datos.
    - `over::AbstractDict`: Diccionario cuyos datos se superpondrán sobre el diccionario base.

    # Devuelve
    - `AbstractDict`: El mismo diccionario `base` ya mutado y actualizado.

    # Excepciones
    - No lanza excepciones propias durante la fusión de datos.
    """

    # Recorre cada par clave-valor del diccionario con los datos a superponer
    for (k, v) in over
        # Si la clave existe en ambos y los dos valores son diccionarios, se fusionan recursivamente
        if v isa AbstractDict && haskey(base, k) && base[k] isa AbstractDict
            deep_merge!(base[k], v)
        else
            # En caso contrario, asigna una copia independiente del valor a la clave correspondiente
            base[k] = deepcopy(v)
        end
    end
    # Devuelve el diccionario base modificado
    return base
end

# Lista que define el orden deseado de las claves en la configuración
const KEY_ORDER = [
    "Duration_hours", "Data_time_step_seconds", "init_epoch",
    "navigation_error_flag", "graphic_display_flag", "image_generation_flag",
    "ignore_fov_for_measurements_flag", "star_csv_flag",
    "camera_config", "photometry_config", "objectives", "observers", "stations_llh",
    "Vmag_lim", "fov_deg", "resolution_x", "resolution_y", "mag0", "flux0", "fwhm_deg",
    "background_level", "sigma_read_noise", "upper_percentile", "lower_percentage", "los_sigma_arcsec",
    "phi_sun", "glare_angle_deg", "Vref_phong", "Fref_phong",
    "SMA", "ECC", "INC", "RAAN", "AOP", "TA",
    "w_deg_s", "w_dir", "obs_dir", "position_sigma_km", "velocity_sigma_km_s"
]

# Convierte la lista anterior en una tabla de consulta rápida. 
# Permite saber inmediatamente qué clave va antes que otra sin tener que recorrer la lista entera.
const KEY_PRIORITY = Dict(k => i for (i, k) in enumerate(KEY_ORDER))

# Claves del esquema en su orden; el resto en orden natural (cam2 antes que cam10).
function sortkey(k::AbstractString)
    """
    Transforma la clave `k` en una tupla de 3 valores para permitir una ordenación personalizable y natural.

    Esta función permite que:
    1. Las claves de alta prioridad (definidas en `KEY_PRIORITY`) aparezcan primero en el JSON.
    2. Las claves con números al final (ej. "cam2" y "cam10") se ordenen por valor numérico real y no alfabéticamente.

    # Argumentos
    - `k::AbstractString`: Clave o texto a clasificar (ej. `"type"`, `"cam2"`, `"cam10"`).

    # Devuelve
    - `Tuple{Int, String, Int}`: Tupla `(prioridad, texto_base, numero_final)`.
      - Elemento 1: Prioridad de la clave (menor número = va más arriba).
      - Elemento 2: Texto antes de los dígitos (ej. `"cam"`).
      - Elemento 3: Valor numérico al final como entero (ej. `2`), o `0` si no hay dígitos.

    # Excepciones
    - `ArgumentError`: Se lanza si falla la conversión del sufijo a número entero (`parse(Int, ...)`).

    # Ejemplo de uso con `sort`
    ```julia
    claves = ["cam10", "type", "cam2"]
    
    # Al aplicar sortkey:
    # "type"  -> (1, "type", 0)           <- Prioridad 1
    # "cam2"  -> (9223372036854775807, "cam", 2)
    # "cam10" -> (9223372036854775807, "cam", 10)
    
    sort(claves, by = sortkey)
    # Resultado ordenado: ["type", "cam2", "cam10"]
    ```
    """

    # 1. Asigna prioridad baja (1, 2...) si está en KEY_PRIORITY, o la máxima posible para enviarla al final
    pri = get(KEY_PRIORITY, k, typemax(Int))
    
    # 2. Separa el texto de los dígitos finales; m es nada si la clave no acaba en números
    m = match(r"^(.*?)(\d+)$", k)
    
    # 3. Devuelve la tupla resultante: si m === nothing no hay número (usa 0), de lo contrario parsea el número como entero
    return m === nothing ? (pri, String(k), 0) : (pri, String(m.captures[1]), parse(Int, m.captures[2]))
end

function write_scalar(io::IO, x)
    """
    Escribe un valor escalar (número, texto, booleano, etc.) formateado en JSON sobre un flujo de entrada/salida (input/output stream, `io`).

    # Argumentos
    - `io::IO`: Canal o flujo de salida (un archivo abierto, la consola `stdout`, un `IOBuffer`, etc.).
    - `x`: Valor escalar que se va a convertir y escribir en formato JSON.

    # Devuelve
    - `Nothing`: La función no retorna ningún valor.

    # Excepciones
    - `ErrorException`: Se lanza si `x` es un número flotante no finito (`NaN` o `Inf`), ya que JSON no los soporta.

    # Ejemplo de uso
    ```julia
    # Escribir un número flotante en la consola
    write_scalar(stdout, 3.1416) # Muestra: 3.1416

    # Escribir una cadena en un archivo de texto
    open("salida.json", "w") do archivo
        write_scalar(archivo, "Hola mundo")
    end
    ```
    """

    # Lanza un error si el valor es NaN o Inf, ya que el estándar JSON no admite estos flotantes
    x isa AbstractFloat && !isfinite(x) && error("Non-finite value cannot be written to JSON.")
    
    # Convierte el escalar a sintaxis JSON y lo escribe directamente en el canal o flujo io
    print(io, JSON3.write(x))
end

function write_value(io::IO, x, level::Int)
    """
    Serializa recursivamente cualquier estructura de datos `x` en formato JSON con sangría formateada (*pretty-print*).

    Organiza diccionarios con sangría de 4 espacios según el nivel de anidamiento (`level`) y aplica `sortkey` 
    para ordenar sus claves. Escribe los vectores en una sola línea compacta y delega los valores escalares a `write_scalar`.

    # Argumentos
    - `io::IO`: Canal o flujo de salida donde se imprimirá la estructura JSON.
    - `x`: Estructura de datos a procesar (diccionario, vector o escalar).
    - `level::Int`: Nivel actual de anidamiento (profundidad) utilizado para aplicar el número de espacios de sangría.

    # Devuelve
    - `Nothing`: La función no retorna ningún valor.

    # Excepciones
    - `ErrorException`: Se lanza indirectamente a través de `write_scalar` si se encuentra un valor flotante `NaN` o `Inf`.

    # Ejemplo de uso
    ```julia
    datos = Dict("type" => "LEO", "cam1" => [1, 2, 3])

    # Imprimir en consola con sangría desde el nivel base (0)
    write_value(stdout, datos, 0)
    # Salida:
    # {
    #     "type": "LEO",
    #     "cam1": [1, 2, 3]
    # }
    ```
    """

    if x isa AbstractDict
        # Si el diccionario no tiene elementos, escribe las llaves vacías {} y finaliza
        isempty(x) && return print(io, "{}")
        
        # Obtiene todas las claves del diccionario y las ordena usando la función personalizada sortkey
        ks = sort!(String.(collect(keys(x))); by = sortkey)
        println(io, "{")
        
        # Recorre cada clave del diccionario aplicando la sangría correspondiente al nivel de profundidad
        for (j, k) in enumerate(ks)
            print(io, "    "^(level + 1), JSON3.write(k), ": ")
            write_value(io, x[k], level + 1)
            println(io, j < length(ks) ? "," : "") # Escribe una coma al final del par salvo si es el último
        end
        print(io, "    "^level, "}")
        
    elseif x isa AbstractVector
        # Para vectores, abre corchetes y escribe todos sus elementos en una misma línea separados por comas
        print(io, "[")
        for (j, v) in enumerate(x)
            j > 1 && print(io, ", ")
            write_scalar(io, v)
        end
        print(io, "]")
        
    else
        # Para valores individuales no compuestos, delega la escritura a la función write_scalar
        write_scalar(io, x)
    end
end

# Guarda el diccionario de configuración en un archivo en disco
function write_config(path::AbstractString, cfg::AbstractDict)
    """
    Escribe un diccionario de configuración en un archivo JSON en disco con formato estructurado.

    Crea automáticamente la estructura de directorios necesaria si no existe previamente en la ruta.

    # Argumentos
    - `path::AbstractString`: Ruta de destino (relativa o absoluta) donde se guardará el archivo JSON.
    - `cfg::AbstractDict`: Diccionario que contiene la configuración a serializar.

    # Devuelve
    - `String`: Devuelve la ruta del archivo escrito tal como se recibió en `path`.

    # Excepciones
    - `SystemError`: Se lanza si ocurren problemas de permisos o de E/S al crear carpetas o escribir el archivo.
    - `ErrorException`: Se lanza (vía `write_value`) si el diccionario contiene valores numéricos no finitos (`NaN` o `Inf`).

    # Ejemplo de uso
    ```julia
    config = Dict("type" => "LEO", "cam1" => Dict("fps" => 30))

    # Guarda el diccionario creando automáticamente las carpetas necesarias
    write_config("configuraciones/simulacion.json", config)
    ```
    """

    # Crea en el sistema de archivos la ruta completa de carpetas contenedoras si no existen
    mkpath(dirname(abspath(path)))
    
    # Abre o crea el archivo para escritura y escribe el contenido con formato JSON
    open(path, "w") do io
        write_value(io, cfg, 0)
        println(io) # Agrega un salto de línea final al terminar de escribir la estructura
    end
    
    # Devuelve la ruta especificada
    return path
end

# -----------------------------------------------------------------------------
# Validación
# -----------------------------------------------------------------------------

function validate_elements(name::AbstractString, el::AbstractDict, ctx)
    """
    Valida la presencia de parámetros, los rangos físicos y la viabilidad orbital de un conjunto de elementos keplerianos.

    # Argumentos
    - `name::AbstractString`: Identificador o nombre descriptivo de la órbita (para mensajes de error).
    - `el::AbstractDict`: Diccionario con los elementos orbitales (debe contener claves como `"SMA"`, `"ECC"`, `"INC"`).
    - `ctx`: Objeto de contexto que proporciona constantes físicas necesarias (ej. `ctx.R` para el radio terrestre).

    # Devuelve
    - `Nothing`: Devuelve `nothing` si todas las validaciones son superadas con éxito.

    # Excepciones
    - `ErrorException`: Se lanza si falta alguna clave obligatoria en `el`.
    - `ErrorException`: Se lanza si el semieje mayor (`SMA`) no es positivo.
    - `ErrorException`: Se lanza si la excentricidad (`ECC`) no está en el rango `[0, 1)`.
    - `ErrorException`: Se lanza si la inclinación (`INC`) no está en `[-180, 180]` grados.
    - `ErrorException`: Se lanza si el perigeo calculado sitúa la órbita por debajo de la superficie terrestre.

    # Ejemplo de uso
    ```julia
    # Contexto simulado con el radio de la Tierra en km
    ctx = (R = 6378.137,)

    elementos = Dict("SMA" => 7000.0, "ECC" => 0.001, "INC" => 98.2, "RAAN" => 0.0, "AOP" => 0.0, "TA" => 0.0)

    # Valida la integridad física de la órbita
    validate_elements("Orbita_SSO", elementos, ctx)
    ```
    """

    # Comprueba que el diccionario posea todas las claves requeridas por la constante ELEMENT_KEYS
    for k in ELEMENT_KEYS
        haskey(el, k) || error("'$name': element '$k' missing (the orbit plugin must return all of $(ELEMENT_KEYS)).")
    end

    # Extrae las variables del semieje mayor, excentricidad e inclinación
    a, e, i = el["SMA"], el["ECC"], el["INC"]

    # Valida que el semieje mayor sea estricatmente positivo
    a > 0           || error("'$name': SMA must be positive (got $a).")

    # Valida que la excentricidad sea apropiada para una órbita elíptica o circular restringida
    0 <= e < 1      || error("'$name': ECC must be in [0,1) (got $e).")

    # Valida que el ángulo de inclinación orbital se encuentre en el rango permitido
    -180 <= i <= 180 || error("'$name': INC must be in [-180,180] deg (got $i).")

    # Calcula la distancia desde el centro de la Tierra al punto más cercano de la órbita (perigeo)
    rp = a * (1 - e)
    
    # Comprueba que el perigeo esté por encima de la superficie terrestre
    rp > ctx.R || error("'$name': perigee radius $(round(rp, digits=1)) km is below the Earth radius.")
    
    # Emite una advertencia si la altitud en el perigeo es menor a 150 km debido a la degradación atmosférica
    rp < ctx.R + 150 && @warn "'$name': perigee altitude < 150 km, the orbit decays very fast."
    
    # Finaliza la función devolviendo nada si la validación es correcta
    return nothing
end

# -----------------------------------------------------------------------------
# Construcción
# -----------------------------------------------------------------------------

function resolve_pointing(group::AbstractDict, member::AbstractDict, name::AbstractString, ctx)
    """
    Resuelve y obtiene el vector o configuración de apuntamiento para un miembro observador.

    Busca y ejecuta la función de apuntamiento registrada en `POINTINGS` según el tipo especificado en el grupo.

    # Argumentos
    - `group::AbstractDict`: Diccionario de configuración del grupo observador (debe contener la clave `"pointing"`).
    - `member::AbstractDict`: Datos o elementos orbitales individuales del miembro del grupo.
    - `name::AbstractString`: Nombre identificador del miembro observador (para contexto de error).
    - `ctx`: Objeto de contexto de simulación pasado a la función de apuntamiento.

    # Devuelve
    - Resultado devuelto por el plugin de apuntamiento invocado (ej. vector o dirección).

    # Excepciones
    - `ErrorException`: Se lanza si `group` no contiene el bloque `"pointing"`.
    - `ErrorException`: Se lanza si el bloque `"pointing"` no define la clave `"type"`.
    - `ErrorException`: Se lanza (vía `lookup`) si el tipo especificado no existe en `POINTINGS`.

    # Ejemplo de uso
    ```julia
    grupo = Dict("pointing" => Dict("type" => "nadir"))
    miembro = Dict("SMA" => 7000.0)
    
    dir_apuntamiento = resolve_pointing(grupo, miembro, "cam1", ctx)
    ```
    """

    # Verifica que el grupo posea el bloque de configuración "pointing"; lanza error si no
    haskey(group, "pointing") || error("Observer group for '$name' needs a 'pointing' block.")
    
    # Extrae el sub-diccionario de configuración del apuntamiento
    pspec = group["pointing"]
    
    # Comprueba que la especificación de apuntamiento incluya la clave "type"
    haskey(pspec, "type") || error("'pointing' block without 'type'.")
    
    # Recupera del registro POINTINGS la función constructora según el tipo y la invoca inmediatamente
    return lookup(POINTINGS, "pointing", String(pspec["type"]))(pspec, member, ctx)
end

function build_section(groups::AbstractDict, ctx, role::Symbol)
    """
    Construye y valida una sección completa de objetos orbitales (observadores o un grupo distinto) a partir de patrones.

    Procesa cada grupo generando sus miembros, asignando nombres únicos, validando sus parámetros
    orbitales y adjuntando las propiedades específicas según su rol (`:observer` u otro).

    # Argumentos
    - `groups::AbstractDict`: Diccionario con la definición de todos los grupos a procesar.
    - `ctx`: Objeto de contexto de simulación utilizado en patrones, validaciones y apuntamiento.
    - `role::Symbol`: Rol de los miembros a construir (`:observer` para observadores o cualquier otro símbolo).

    # Devuelve
    - `Dict{String,Any}`: Diccionario procesado con el nombre de cada miembro como clave y sus parámetros construidos como valor.

    # Excepciones
    - `ErrorException`: Se lanza si el patrón invocado no genera miembros (`members` está vacío).
    - `ErrorException`: Se lanza si se genera un nombre de miembro duplicado en el diccionario de salida.
    - `ErrorException`: Se lanza (vía `validate_elements`) si algún miembro no cumple con las restricciones orbitales físicas.
    - `ErrorException`: Se lanza (vía `resolve_pointing`) si falla la resolución de apuntamiento cuando `role === :observer`.

    # Ejemplo de uso
    ```julia
    grupos = Dict(
        "cam" => Dict(
            "pattern" => "single",
            "pointing" => Dict("type" => "nadir"),
            "SMA" => 7000.0, "ECC" => 0.001, "INC" => 98.0,
            "RAAN" => 0.0, "AOP" => 0.0, "TA" => 0.0
        )
    )

    seccion_obs = build_section(grupos, ctx, :observer)
    ```
    """

    # Inicializa el diccionario de salida para almacenar los miembros procesados
    out = Dict{String,Any}()

    # Recorre cada grupo de configuración especificado en la entrada
    for (gname, group) in groups
        # Obtiene la clave del patrón de distribución (usa "single" por defecto si no existe)
        pattern = String(get(group, "pattern", "single"))
        
        # Genera la lista de miembros llamando a la función del patrón registrada en PATTERNS
        members = lookup(PATTERNS, "pattern", pattern)(group, ctx)
        
        # Lanza un error si la generación del patrón no produjo ningún miembro
        isempty(members) && error("Group '$gname': pattern '$pattern' returned no members.")

        # Guarda la cantidad total de miembros generados para controlar el nombrado
        n = length(members)

        # Procesa cada miembro individual generado por el patrón
        for (idx, el) in enumerate(members)
            # Determina el nombre: si hay solo 1 usa el nombre del grupo; si hay varios añade un índice (ej. cam1, cam2)
            name = n == 1 ? String(gname) : string(gname, idx)
            
            # Asegura que el nombre no esté repetido en el resultado final
            haskey(out, name) && error("Duplicate name '$name'.")

            # Valida los parámetros orbitales keplerianos del miembro
            validate_elements(name, el, ctx)

            # Crea el diccionario base del miembro copiando únicamente las claves requeridas
            entry = Dict{String,Any}(k => el[k] for k in ELEMENT_KEYS)

            # Agrega los campos específicos dependiendo de si el rol es :observer o de otro tipo
            if role === :observer
                # Asigna la dirección de apuntamiento y las desviaciones típicas de posición y velocidad
                entry["obs_dir"] = resolve_pointing(group, el, name, ctx)
                entry["position_sigma_km"]   = get(group, "position_sigma_km",   ctx.member_fields["position_sigma_km"])
                entry["velocity_sigma_km_s"] = get(group, "velocity_sigma_km_s", ctx.member_fields["velocity_sigma_km_s"])
            else
                # Asigna velocidad angular y su dirección para roles que no sean observadores
                entry["w_deg_s"] = get(group, "w_deg_s", ctx.member_fields["w_deg_s"])
                entry["w_dir"]   = get(group, "w_dir",   ctx.member_fields["w_dir"])
            end

            # Registra la entrada construida en el diccionario final usando su nombre único
            out[name] = entry
        end
    end

    # Devuelve el diccionario final compilado con todos los miembros procesados
    return out
end

function generate_config(spec::AbstractDict; spec_dir::AbstractString = pwd())
    """
    Genera un diccionario de configuración completo y validado a partir de una especificación dada.

    Fusiona la plantilla por defecto con plantillas externas y el escenario indicado, construye las secciones
    de objetivos y observadores utilizando sus respectivos patrones y valida la consistencia global.

    # Argumentos
    - `spec::AbstractDict`: Diccionario con la especificación requerida (puede incluir `"template"`, `"scenario"`, `"objectives"`, etc.).

    # Argumentos de palabra clave (kwargs)
    - `spec_dir::AbstractString = pwd()`: Directorio base para resolver rutas relativas a archivos de plantilla.

    # Devuelve
    - `Dict{String,Any}`: Diccionario con la configuración completa compilada y validada.

    # Excepciones
    - `ErrorException`: Se lanza si no se define al menos un objetivo en `spec`.
    - `ErrorException`: Se lanza si no se define al menos un observador en `spec`.
    - `ErrorException`: Se lanza (vía `read_json` o `build_section`) si fallan la lectura de archivos o la validación física orbital.
    ```
    """

    # Copia la plantilla base predeterminada y la estructura de campos de miembros para no mutar los originales
    cfg = deepcopy(DEFAULT_TEMPLATE)
    member_fields = deepcopy(DEFAULT_MEMBER_FIELDS)

    # Si se especifica una plantilla JSON externa, la lee y la fusiona con la configuración base
    if haskey(spec, "template")
        tpath = String(spec["template"])
        isabspath(tpath) || (tpath = normpath(joinpath(spec_dir, tpath)))
        deep_merge!(cfg, read_json(tpath))
        delete!(cfg, "objectives")
        delete!(cfg, "observers")
    end

    # Superpone la configuración específica del escenario y los campos de miembros si están presentes
    haskey(spec, "scenario") && deep_merge!(cfg, spec["scenario"])
    haskey(spec, "member_fields") && deep_merge!(member_fields, spec["member_fields"])

    # Crea la tupla de contexto con fecha de época, constantes orbitales y campos asignados
    ctx = (epoch = DateTime(String(cfg["init_epoch"])), mu = μ, R = R, member_fields = member_fields)

    # Construye la sección de objetivos y lanza un error si la lista resultante está vacía
    objectives = build_section(get(spec, "objectives", Dict{String,Any}()), ctx, :objective)
    isempty(objectives) && error("The spec needs at least one objective.")

    # Añade los objetivos al contexto antes de procesar observadores para permitir apuntamientos referenciados a metas
    ctx_obs = merge(ctx, (objectives = objectives,))
    observers = build_section(get(spec, "observers", Dict{String,Any}()), ctx_obs, :observer)
    isempty(observers) && error("The spec needs at least one observer.")

    # Asigna las estructuras de objetivos y observadores procesadas al diccionario final
    cfg["objectives"] = objectives
    cfg["observers"]  = observers

    # Devuelve el diccionario de configuración completamente generado
    return cfg
end

function generate_config_file(spec_path::AbstractString, out_path::AbstractString)
    """
    Lee un archivo JSON de especificación, procesa la configuración y la guarda en un nuevo archivo JSON en disco.

    # Argumentos
    - `spec_path::AbstractString`: Ruta del archivo JSON de entrada con las especificaciones.
    - `out_path::AbstractString`: Ruta de salida donde se guardará el archivo JSON resultante.

    # Devuelve
    - `String`: Devuelve la ruta absoluta o relativa del archivo generado especificada en `out_path`.

    # Excepciones
    - `ErrorException`: Se lanza (vía `read_json` o `generate_config`) si el archivo no existe, no es un JSON válido o la configuración no cumple los requisitos.
    - `SystemError`: Se lanza si ocurren fallos al escribir el archivo de salida en disco.    
    """

    # Obtiene la ruta absoluta del archivo de especificación para resolver rutas relativas asociadas
    spec_path = abspath(spec_path)
    
    # Lee e interpreta el archivo JSON de especificación
    spec = read_json(spec_path)
    
    # Genera el diccionario de configuración tomando como directorio base la ubicación de spec_path
    cfg = generate_config(spec; spec_dir = dirname(spec_path))
    
    # Escribe la configuración generada en la ruta de salida indicada
    write_config(out_path, cfg)

    # Imprime en la consola un resumen de la generación indicando las claves ordenadas
    println("Config generated: $out_path")
    println("  objectives: ", join(sort!(collect(keys(cfg["objectives"])); by = sortkey), ", "))
    println("  observers:  ", join(sort!(collect(keys(cfg["observers"]));  by = sortkey), ", "))
    
    # Devuelve la ruta del archivo generado
    return out_path
end
