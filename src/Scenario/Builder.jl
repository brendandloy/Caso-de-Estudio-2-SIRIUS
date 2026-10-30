# =============================================================================
# Builder.jl
#
# Núcleo del generador: convierte una spec en el JSON completo de escenario.
# Se incluye desde Scenario.jl, por lo que comparte su espacio de nombres
# (ve getreq, getopt, lookup, los registros, etc.).
#
# Constantes:
#   ELEMENT_KEYS, DEFAULT_TEMPLATE, DEFAULT_MEMBER_FIELDS, SPEC_KEYS, FREE_KEY_BLOCKS,
#   KEY_ORDER, KEY_PRIORITY
#
# Si se añade o cambia un valor de DEFAULT_TEMPLATE o DEFAULT_MEMBER_FIELDS, hay que
# actualizar también config/specs/plantilla.json, que los lista todos.
#
# Funciones:
#   read_json, deep_merge!          : lectura de JSON y fusión recursiva de diccionarios
#   sortkey, write_scalar,
#   write_value, write_config       : escritura del JSON con orden de claves estable
#   validate_elements               : comprobaciones físicas de cada miembro
#   check_known_keys                : rechaza claves desconocidas en scenario / member_fields
#   resolve_pointing                : obtiene obs_dir con el plugin de apuntado
#   resolve_attitude                : bloque "attitude" de cada objetivo (validado con el módulo Attitude)
#   build_section                  : construye objectives u observers
#   generate_config                 : spec -> Dict completo de escenario
#   generate_config_file            : spec.json -> JSON en disco
#
# Prioridad de valores (de menor a mayor):
#   - Escenario: DEFAULT_TEMPLATE < archivo "template" < bloque "scenario".
#   - Miembros:  DEFAULT_MEMBER_FIELDS < bloque "member_fields" < campo del grupo.
#   - Actitud:   sin valor por defecto; cada grupo de objetivos la define con su
#                bloque "attitude" (las claves antiguas w_deg_s / w_dir dan error).
# =============================================================================

# Los seis elementos keplerianos que lee el simulador
const ELEMENT_KEYS = ("SMA", "ECC", "INC", "RAAN", "AOP", "TA")

# Valores de fábrica del escenario (los de tus config_*.json); se pisan con "template" y "scenario"
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

# Valores de fábrica de los campos propios de cada miembro (cámaras)
# La actitud de los objetivos no tiene valor de fábrica: es obligatoria en cada grupo
const DEFAULT_MEMBER_FIELDS = Dict{String,Any}(
    "position_sigma_km"   => [0.015, 0.015, 0.015],
    "velocity_sigma_km_s" => [0.00002, 0.00002, 0.00002]
)

# Claves admitidas en el nivel superior de una spec
const SPEC_KEYS = ("template", "seed", "scenario", "member_fields", "objectives", "observers")

# Bloques de DEFAULT_TEMPLATE cuyas claves son nombres libres (no se comprueban)
const FREE_KEY_BLOCKS = ("stations_llh",)

# -----------------------------------------------------------------------------
# JSON: lectura, fusión y escritura (con orden de claves estable)
# -----------------------------------------------------------------------------

"""
    read_json(path) -> Dict{String,Any}

Lee un archivo JSON y lo devuelve como diccionario.

# Argumentos
- `path::AbstractString`: ruta del archivo JSON.

# Devuelve
- `Dict{String,Any}` con el contenido del archivo.

# Excepciones
- `ErrorException`: si el archivo no existe.
- `ErrorException`: si el JSON no se puede interpretar (p. ej. llaves sobrantes al final).

# Ejemplo
```julia
spec = read_json("config/specs/ejemplo_leo.json")
keys(spec)        # "scenario", "objectives", "observers"
```
"""
function read_json(path::AbstractString)
    isfile(path) || error("File not found: $path")
    try
        return JSON3.read(read(path, String), Dict{String,Any})
    catch err
        error("Could not parse '$path' (revisa llaves sobrantes al final): " * sprint(showerror, err))
    end
end

"""
    deep_merge!(base, over) -> base

Fusiona `over` dentro de `base` de forma recursiva: los diccionarios se
mezclan clave a clave y cualquier otro valor de `over` sustituye al de `base`.

Por qué existe: permite cambiar un solo valor anidado (p. ej. `fov_deg`) sin
repetir todo el bloque `camera_config` en la spec.

# Argumentos
- `base::AbstractDict`: diccionario que se modifica.
- `over::AbstractDict`: diccionario cuyos valores tienen prioridad.

# Devuelve
- `base`, ya modificado.

# Excepciones
- Ninguna.

# Ejemplo
```julia
base = Dict("camera_config" => Dict("fov_deg" => 27, "mag0" => 8.0))
over = Dict("camera_config" => Dict("fov_deg" => 20))
deep_merge!(base, over)     # camera_config = Dict("fov_deg" => 20, "mag0" => 8.0)
```
"""
function deep_merge!(base::AbstractDict, over::AbstractDict)
    for (k, v) in over
        if v isa AbstractDict && haskey(base, k) && base[k] isa AbstractDict
            deep_merge!(base[k], v)         # ambos son diccionarios: se mezclan
        else
            base[k] = deepcopy(v)           # se copia para no compartir memoria con `over`
        end
    end
    return base
end

# Orden de claves del JSON de salida (el simulador no lo exige, pero queda legible)
const KEY_ORDER = [
    "Duration_hours", "Data_time_step_seconds", "init_epoch",
    "navigation_error_flag", "graphic_display_flag", "image_generation_flag",
    "ignore_fov_for_measurements_flag", "star_csv_flag",
    "camera_config", "photometry_config", "objectives", "observers", "stations_llh",
    "Vmag_lim", "fov_deg", "resolution_x", "resolution_y", "mag0", "flux0", "fwhm_deg",
    "background_level", "sigma_read_noise", "upper_percentile", "lower_percentage", "los_sigma_arcsec",
    "phi_sun", "glare_angle_deg", "Vref_phong", "Fref_phong",
    "SMA", "ECC", "INC", "RAAN", "AOP", "TA",
    "attitude", "obs_dir", "position_sigma_km", "velocity_sigma_km_s",
    # Dentro del bloque "attitude": primero el modo y luego sus parámetros
    "type", "w_deg_s", "w_dir", "seed", "frame", "q", "euler_deg", "sequence",
    "primary", "secondary", "body_axis", "target"
]

# Clave -> posición en KEY_ORDER
const KEY_PRIORITY = Dict(k => i for (i, k) in enumerate(KEY_ORDER))

"""
    sortkey(k) -> Tuple{Int,String,Int}

Clave de ordenación de un nombre: primero las claves del esquema (en el
orden de KEY_ORDER) y luego el resto en orden natural (cam2 antes que cam10).

Por qué existe: Julia no garantiza el orden de las claves de un `Dict`, y sin
esto el JSON saldría desordenado (o con cam10 antes que cam2).

# Argumentos
- `k::AbstractString`: nombre de la clave.

# Devuelve
- `Tuple{Int,String,Int}`: (prioridad, texto sin dígitos finales, número final).

# Excepciones
- Ninguna.

# Ejemplo
```julia
sort(["cam10", "cam2", "SMA"]; by = sortkey)     # ["SMA", "cam2", "cam10"]
```
"""
function sortkey(k::AbstractString)
    pri = get(KEY_PRIORITY, k, typemax(Int))      # las claves del esquema van antes
    m = match(r"^(.*?)(\d+)$", k)                 # separa el texto del número final
    return m === nothing ? (pri, String(k), 0) : (pri, String(m.captures[1]), parse(Int, m.captures[2]))
end

"""
    write_scalar(io, x)

Escribe un valor escalar como JSON en `io`.

Por qué existe: JSON no admite NaN ni Inf; esta función lo detecta antes de
escribir un JSON que el simulador no podría leer.

# Argumentos
- `io::IO`: flujo de salida.
- `x`: número, texto o booleano.

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si `x` es un decimal no finito (NaN o Inf).

# Ejemplo
```julia
write_scalar(stdout, 2.5)        # imprime 2.5
write_scalar(stdout, NaN)        # error: Non-finite value cannot be written to JSON.
```
"""
function write_scalar(io::IO, x)
    x isa AbstractFloat && !isfinite(x) && error("Non-finite value cannot be written to JSON.")
    print(io, JSON3.write(x))
end

"""
    write_value(io, x, level)

Escribe recursivamente un valor como JSON con sangría de 4 espacios por nivel.
Los diccionarios salen con las claves ordenadas por `sortkey` y los vectores
en una sola línea.

Por qué existe: `JSON3.pretty` pone cada número de un vector en su propia
línea y no ordena las claves; este formato imita el de tus JSON escritos a mano.

# Argumentos
- `io::IO`: flujo de salida.
- `x`: diccionario, vector o escalar.
- `level::Int`: nivel de anidamiento actual (controla la sangría).

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si hay un decimal no finito en `x`.

# Ejemplo
```julia
write_value(stdout, Dict("b" => [1, 2], "a" => 1.5), 0)
# {
#     "a": 1.5,
#     "b": [1, 2]
# }
```
"""
function write_value(io::IO, x, level::Int)
    if x isa AbstractDict
        isempty(x) && return print(io, "{}")
        ks = sort!(String.(collect(keys(x))); by = sortkey)
        println(io, "{")
        for (j, k) in enumerate(ks)
            print(io, "    "^(level + 1), JSON3.write(k), ": ")
            write_value(io, x[k], level + 1)
            println(io, j < length(ks) ? "," : "")      # sin coma tras la última clave
        end
        print(io, "    "^level, "}")
    elseif x isa AbstractVector
        print(io, "[")
        for (j, v) in enumerate(x)
            j > 1 && print(io, ", ")
            write_scalar(io, v)
        end
        print(io, "]")
    else
        write_scalar(io, x)
    end
end

"""
    write_config(path, cfg) -> String

Escribe el escenario completo como JSON en disco, creando las carpetas necesarias.

# Argumentos
- `path::AbstractString`: ruta del archivo de salida.
- `cfg::AbstractDict`: escenario completo.

# Devuelve
- `String`: la ruta escrita.

# Excepciones
- `ErrorException`: si `cfg` contiene un decimal no finito.
- `SystemError`: si no se puede crear la carpeta o escribir el archivo.

# Ejemplo
```julia
cfg = generate_config(spec)
write_config("config/generated/mi_escenario.json", cfg)
```
"""
function write_config(path::AbstractString, cfg::AbstractDict)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        write_value(io, cfg, 0)
        println(io)
    end
    return path
end

# -----------------------------------------------------------------------------
# Validación
# -----------------------------------------------------------------------------

"""
    validate_elements(name, el, ctx)

Comprueba que los elementos de un miembro son físicamente razonables.

# Argumentos
- `name::AbstractString`: nombre del miembro, para los mensajes de error.
- `el::AbstractDict`: elementos devueltos por el plugin de órbita.
- `ctx`: contexto del Builder (usa `ctx.R`).

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si falta algún elemento, `SMA <= 0`, `ECC` fuera de [0,1),
  `INC` fuera de [-180,180] o el perigeo queda bajo la superficie terrestre.
- Emite un aviso (no un error) si la altitud del perigeo es menor de 150 km.

# Ejemplo
```julia
el = Dict("SMA" => 7000.0, "ECC" => 0.001, "INC" => 90.0, "RAAN" => 0.0, "AOP" => 0.0, "TA" => 0.0)
validate_elements("cam1", el, ctx)        # correcto, devuelve nothing
el["SMA"] = 6000.0
validate_elements("cam1", el, ctx)        # error: perigee radius ... is below the Earth radius
```
"""
function validate_elements(name::AbstractString, el::AbstractDict, ctx)
    for k in ELEMENT_KEYS
        haskey(el, k) || error("'$name': element '$k' missing (the orbit plugin must return all of $(ELEMENT_KEYS)).")
    end

    a, e, i = el["SMA"], el["ECC"], el["INC"]
    a > 0            || error("'$name': SMA must be positive (got $a).")
    0 <= e < 1       || error("'$name': ECC must be in [0,1) (got $e).")
    -180 <= i <= 180 || error("'$name': INC must be in [-180,180] deg (got $i).")

    rp = a * (1 - e)                                # radio del perigeo
    rp > ctx.R || error("'$name': perigee radius $(round(rp, digits=1)) km is below the Earth radius.")
    rp < ctx.R + 150 && @warn "'$name': perigee altitude < 150 km, the orbit decays very fast."
    return nothing
end

# -----------------------------------------------------------------------------
# Construcción
# -----------------------------------------------------------------------------

"""
    resolve_pointing(group, member, name, ctx) -> Vector{Float64}

Obtiene la dirección de apuntado (`obs_dir`) de un observador usando el
plugin de apuntado indicado en el bloque `pointing` de su grupo.

# Argumentos
- `group::AbstractDict`: grupo de la spec al que pertenece el observador.
- `member::AbstractDict`: elementos orbitales del observador.
- `name::AbstractString`: nombre del observador, para los mensajes de error.
- `ctx`: contexto del Builder.

# Devuelve
- `Vector{Float64}` de 3 componentes.

# Excepciones
- `ErrorException`: si el grupo no tiene `pointing`, o `pointing` no tiene `type`.
- `ErrorException`: si el tipo de apuntado no está registrado.

# Ejemplo
```julia
group = Dict("pointing" => Dict("type" => "fixed", "dir" => [0, 2, 0]))
resolve_pointing(group, member, "cam1", ctx)      # [0.0, 1.0, 0.0]
```
"""
function resolve_pointing(group::AbstractDict, member::AbstractDict, name::AbstractString, ctx)
    haskey(group, "pointing") || error("Observer group for '$name' needs a 'pointing' block.")
    pspec = group["pointing"]
    haskey(pspec, "type") || error("'pointing' block without 'type'.")
    return lookup(POINTINGS, "pointing", String(pspec["type"]))(pspec, member, ctx)
end

"""
    resolve_attitude(group, name, ctx) -> Dict{String,Any}

Obtiene el bloque `attitude` de un objetivo a partir de su grupo de la spec y
lo valida con el mismo código que usa el simulador (`Attitude.parse_attitude`).

Reglas:
- El grupo debe tener un bloque `attitude` (mismas reglas que el simulador,
  ver `Attitude.attitude_block`). Las claves antiguas `w_deg_s` / `w_dir` dan error.
- Modo `random`: cada miembro recibe su propia semilla, sacada de un flujo
  aleatorio de la spec separado del de las órbitas (cambiar la actitud no altera
  las órbitas aleatorias). Por eso la spec no admite `seed` dentro del bloque.

# Argumentos
- `group::AbstractDict`: grupo de la spec al que pertenece el objetivo.
- `name::AbstractString`: nombre del objetivo, para los mensajes de error.
- `ctx`: contexto del Builder (usa `ctx.rng_attitude`).

# Devuelve
- `Dict{String,Any}`: bloque `attitude` listo para el JSON (copia; no modifica la spec).

# Excepciones
- `ErrorException`: si el grupo no tiene bloque `attitude`, usa las claves antiguas
  `w_deg_s` / `w_dir` o pone `seed` en un bloque `random`.
- Cualquier error de validación de `Attitude.parse_attitude`.

# Ejemplo
```julia
resolve_attitude(Dict("attitude" => Dict("type" => "spin", "w_deg_s" => 0)), "deb", ctx)
# Dict("type" => "spin", "w_deg_s" => 0)
resolve_attitude(Dict("attitude" => Dict("type" => "random")), "deb1", ctx)
# Dict("type" => "random", "seed" => 2861592718)
```
"""
function resolve_attitude(group::AbstractDict, name::AbstractString, ctx)
    # Copia del bloque (el mismo grupo da varios miembros); error si falta o usa w_deg_s/w_dir
    block = attitude_block(group, name)

    # random: una semilla distinta por miembro, reproducible con la "seed" de la spec
    if get(block, "type", nothing) == "random"
        haskey(block, "seed") && error(
            "Objective '$name': do not put 'seed' in a 'random' attitude of a spec; " *
            "each member gets its own seed from the spec's top-level 'seed'.")
        block["seed"] = Int(rand(ctx.rng_attitude, UInt32))
    end

    parse_attitude(block, name)       # solo valida; el simulador lo volverá a leer del JSON
    return block
end

"""
    build_section(groups, ctx, role) -> Dict{String,Any}

Construye la sección `objectives` u `observers` del escenario: expande cada
grupo con su patrón, valida los miembros y añade los campos propios del rol.

# Argumentos
- `groups::AbstractDict`: grupos de la spec (nombre del grupo -> definición).
- `ctx`: contexto del Builder (incluye `ctx.member_fields`).
- `role::Symbol`: `:observer` (añade `obs_dir` y sigmas) u otro valor, que se trata como objetivo (añade `attitude`).

# Devuelve
- `Dict{String,Any}`: nombre del miembro -> entrada lista para el JSON.
  Nombres: `grupo_etiqueta` si el patrón da `_label` (p. ej. NORAD de un TLE),
  `grupo` si hay un solo miembro y `grupo1..grupoN` si hay varios.

# Excepciones
- `ErrorException`: si el patrón no existe o no devuelve miembros.
- `ErrorException`: si hay nombres duplicados o un miembro no supera la validación.
- `ErrorException`: si un observador no tiene un bloque `pointing` válido.
- `ErrorException`: si un objetivo no tiene una actitud válida (ver `resolve_attitude`).

# Ejemplo
```julia
groups = Dict("cam" => Dict(
    "orbit"    => Dict("type" => "GEO", "longitude_deg" => 0.0),
    "pointing" => Dict("type" => "fixed", "dir" => [0, 1, 0])))
observers = build_section(groups, ctx, :observer)     # Dict con la clave "cam"
```
"""
function build_section(groups::AbstractDict, ctx, role::Symbol)
    out = Dict{String,Any}()

    # Grupos en orden alfabético: con órbitas aleatorias, la misma semilla da el mismo escenario
    for gname in sort!(collect(keys(groups)))
        group = groups[gname]
        pattern = String(get(group, "pattern", "single"))      # "single" si no se indica
        members = lookup(PATTERNS, "pattern", pattern)(group, ctx)
        isempty(members) && error("Group '$gname': pattern '$pattern' returned no members.")

        n = length(members)

        for (idx, el) in enumerate(members)
            # Si el patrón da una etiqueta (p. ej. NORAD de un TLE): grupo_etiqueta (deb_25544)
            # Si no: un miembro usa el nombre del grupo; varios, nombre + índice (cam1, cam2, ...)
            name = haskey(el, "_label") ? string(gname, "_", el["_label"]) :
                   n == 1               ? String(gname) : string(gname, idx)
            haskey(out, name) && error("Duplicate name '$name'.")

            validate_elements(name, el, ctx)

            entry = Dict{String,Any}(k => el[k] for k in ELEMENT_KEYS)

            if role === :observer
                # Observadores: apuntado y errores de navegación (el campo del grupo gana al valor general)
                entry["obs_dir"] = resolve_pointing(group, el, name, ctx)
                entry["position_sigma_km"]   = get(group, "position_sigma_km",   ctx.member_fields["position_sigma_km"])
                entry["velocity_sigma_km_s"] = get(group, "velocity_sigma_km_s", ctx.member_fields["velocity_sigma_km_s"])
            else
                # Objetivos: bloque de actitud (obligatorio en el grupo)
                entry["attitude"] = resolve_attitude(group, name, ctx)
            end

            out[name] = entry
        end
    end

    return out
end

"""
    check_known_keys(block, ref, who)

Comprueba, de forma recursiva, que todas las claves de `block` existen en el
diccionario de referencia `ref` (p. ej. `DEFAULT_TEMPLATE`).

Por qué existe: una errata en `scenario` o `member_fields` (p. ej.
`"image_generation_flg"`) se copiaría al JSON sin efecto y el simulador usaría
el valor por defecto sin avisar. La lista de claves válidas está en
`config/specs/plantilla.json`.

# Argumentos
- `block::AbstractDict`: bloque a comprobar.
- `ref::AbstractDict`: bloque de referencia con todas las claves admitidas.
- `who::AbstractString`: nombre del bloque, para los mensajes de error.

# Devuelve
- `nothing`.

# Excepciones
- `ErrorException`: si aparece una clave que no está en `ref` (el mensaje lista las admitidas).

# Ejemplo
```julia
check_known_keys(Dict("image_generation_flg" => true), DEFAULT_TEMPLATE, "scenario")
# error: 'scenario': unknown parameter 'image_generation_flg'. Allowed: ...
```
"""
function check_known_keys(block::AbstractDict, ref::AbstractDict, who::AbstractString)
    check_keys(block, keys(ref); who = who)

    # Se baja a los sub-bloques (camera_config, photometry_config...), salvo los de nombres libres
    for (k, v) in block
        kk = String(k)
        (v isa AbstractDict && ref[kk] isa AbstractDict && !(kk in FREE_KEY_BLOCKS)) || continue
        check_known_keys(v, ref[kk], "$who.$kk")
    end
    return nothing
end

"""
    generate_config(spec; spec_dir=pwd()) -> Dict{String,Any}

Convierte una spec en el diccionario completo de escenario. Mezcla los
valores de fábrica, el archivo `template` (opcional) y el bloque `scenario`,
y construye `objectives` y `observers`.

# Argumentos
- `spec::AbstractDict`: spec ya leída (`template`, `scenario`, `member_fields`, `seed`, `objectives`, `observers`).
- `spec_dir::AbstractString`: carpeta de la spec; resuelve rutas relativas de `template`.

# Devuelve
- `Dict{String,Any}` con el escenario completo, en el esquema de Environment.load_config.

# Excepciones
- `ErrorException`: si la spec, `scenario`, el `template` o `member_fields` tienen una clave desconocida.
- `ErrorException`: si el `template` no existe o no se puede leer.
- `ErrorException`: si la spec no tiene al menos un objetivo y un observador.
- `ErrorException`: cualquiera de los errores de `build_section` y de los plugins.

# Ejemplo
```julia
spec = read_json("config/specs/ejemplo_sso_geo.json")
cfg = generate_config(spec)
cfg["observers"]["cam1"]["INC"]        # ≈ 98.19 (SSO a 700 km)
```
"""
function generate_config(spec::AbstractDict; spec_dir::AbstractString = pwd())
    check_keys(spec, SPEC_KEYS; who = "spec")             # p. ej. "scenarios" por "scenario"

    cfg = deepcopy(DEFAULT_TEMPLATE)                      # valores de fábrica del escenario
    member_fields = deepcopy(DEFAULT_MEMBER_FIELDS)       # valores de fábrica de cada miembro

    if haskey(spec, "template")
        tpath = String(spec["template"])
        isabspath(tpath) || (tpath = normpath(joinpath(spec_dir, tpath)))
        deep_merge!(cfg, read_json(tpath))
        # Del template solo interesa lo común; los miembros los define la spec
        delete!(cfg, "objectives")
        delete!(cfg, "observers")
    end

    haskey(spec, "scenario") && deep_merge!(cfg, spec["scenario"])
    haskey(spec, "member_fields") && deep_merge!(member_fields, spec["member_fields"])

    # Erratas en el template, en "scenario" o en "member_fields" darían valores por defecto sin avisar
    check_known_keys(cfg, DEFAULT_TEMPLATE, "scenario")
    check_known_keys(member_fields, DEFAULT_MEMBER_FIELDS, "member_fields")

    # Semilla de los números aleatorios: la de la spec o una nueva, que se imprime para poder repetir
    if haskey(spec, "seed")
        seed = getint(spec, "seed")
        seed >= 0 || error("'seed' must be >= 0 (got $seed).")
    else
        seed = Int(rand(RandomDevice(), UInt32))
        println("No 'seed' in spec; using seed = $seed (add \"seed\": $seed to the spec to repeat this scenario).")
    end

    # spec_dir permite a los patrones resolver rutas relativas a la spec (p. ej. archivos TLE)
    # rng es el generador aleatorio común a todos los plugins (p. ej. random_LEO)
    # rng_attitude da las semillas de las actitudes "random"; va aparte para no alterar las órbitas
    ctx = (epoch = DateTime(String(cfg["init_epoch"])), mu = μ, R = R,
           member_fields = member_fields, spec_dir = spec_dir, rng = MersenneTwister(seed),
           rng_attitude = MersenneTwister(seed + 1))

    objectives = build_section(get(spec, "objectives", Dict{String,Any}()), ctx, :objective)
    isempty(objectives) && error("The spec needs at least one objective.")

    # Los apuntados pueden necesitar conocer los objetivos (p. ej. "at_target")
    ctx_obs = merge(ctx, (objectives = objectives,))
    observers = build_section(get(spec, "observers", Dict{String,Any}()), ctx_obs, :observer)
    isempty(observers) && error("The spec needs at least one observer.")

    cfg["objectives"] = objectives
    cfg["observers"]  = observers

    return cfg
end

"""
    generate_config_file(spec_path, out_path) -> String

Lee una spec desde disco, genera el escenario completo y lo escribe como
JSON. Imprime los nombres de los objetivos y observadores generados.

# Argumentos
- `spec_path::AbstractString`: ruta de la spec.
- `out_path::AbstractString`: ruta del JSON de salida.

# Devuelve
- `String`: la ruta del JSON escrito.

# Excepciones
- `ErrorException`: si la spec no existe, no se puede leer o no es válida.
- `SystemError`: si no se puede escribir el archivo de salida.

# Ejemplo
```julia
generate_config_file("config/specs/ejemplo_sso_geo.json", "config/generated/ejemplo_sso_geo.json")
```
"""
function generate_config_file(spec_path::AbstractString, out_path::AbstractString)
    spec_path = abspath(spec_path)
    spec = read_json(spec_path)
    cfg = generate_config(spec; spec_dir = dirname(spec_path))
    write_config(out_path, cfg)

    println("Config generated: $out_path")
    println("  objectives: ", join(sort!(collect(keys(cfg["objectives"])); by = sortkey), ", "))
    println("  observers:  ", join(sort!(collect(keys(cfg["observers"]));  by = sortkey), ", "))
    return out_path
end