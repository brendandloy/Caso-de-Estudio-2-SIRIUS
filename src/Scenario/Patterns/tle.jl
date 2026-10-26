# =============================================================================
# tle.jl
#
# Patrón TLE: crea un miembro por cada TLE de un archivo, propagado con SGP4
# hasta la init_epoch del escenario. Registra el patrón "tle".
#
#   "deb": {
#     "pattern": "tle",
#     "file":    "../tle/catalogo.tle",
#     "select":  {"norad_ids": [25544, 47699]},
#     "attitude": {"type": "random", "w_range_deg_s": [0, 5]}
#   }
#
# - "file": archivo con TLEs de 2 o 3 líneas (con nombre). Ruta absoluta o relativa
#   a la carpeta de la spec (igual que "template"): "../tle/a.tle" se interpreta desde
#   donde está el JSON de la spec, no desde donde se lanza Julia.
# - "select" (opcional; sin él se usan todos los TLEs del archivo), una sola clave:
#     "norad_ids": [25544, ...]     por número NORAD
#     "names":     ["ISS (ZARYA)"]  por nombre (línea 0 del TLE)
#     "range":     [1, 50]          por posición en el archivo (1-based, ambos incluidos)
# - Nombres de los miembros: grupo_NORAD (deb_25544), los pone el Builder con "_label".
#
# Cadena de conversión de cada TLE:
#   1) SGP4 (WGS84) desde la época del TLE hasta init_epoch -> r, v en TEME (km, km/s)
#   2) TEME -> J2000 (IAU-76/FK5, sin EOP) con la matriz de rotación en init_epoch
#   3) cartesian_to_kepler con el μ del simulador -> elementos osculadores
#
# Elementos osculadores: los de la elipse kepleriana que seguiría el objeto si en ese
# instante cesaran las perturbaciones (J2, arrastre...). Son tangentes a la trayectoria
# real solo en ese instante y cambian con el tiempo; los de un TLE son medios de SGP4.
#
# Por qué el paso 2: TEME es un sistema "de fecha" y en 2023 está girado ~0.27°
# respecto a J2000; sin la rotación el error de posición en LEO es de ~30 km.
#
# Requiere en el Project.toml: SatelliteToolboxTle, SatelliteToolboxSgp4 y
# SatelliteToolboxTransformations.
#
# Funciones: load_tles, select_tles, tle_to_elements, tle_pattern
# =============================================================================

# Importación explícita: evita choques de nombres con Constants (μ, R...)
using SatelliteToolboxTle: TLE, read_tles_from_file
using SatelliteToolboxSgp4: sgp4_init, sgp4!
using SatelliteToolboxTransformations: r_eci_to_eci, TEME, J2000

# Distancia máxima (días) entre la época del TLE y init_epoch antes de avisar
const TLE_MAX_DAYS = 3.0

"""
    load_tles(path) -> Vector{TLE}

Lee todos los TLEs de un archivo (formato de 2 o 3 líneas).

# Argumentos
- `path::AbstractString`: ruta del archivo de TLEs.

# Devuelve
- `Vector{TLE}` en el orden del archivo. Los TLEs inválidos se descartan con un aviso
  (lo hace `read_tles_from_file`).

# Excepciones
- `ErrorException`: si el archivo no existe o no contiene ningún TLE válido.

# Ejemplo
```julia
tles = load_tles("config/tle/ejemplo.tle")
tles[1].satellite_number        # 47699
```
"""
function load_tles(path::AbstractString)
    isfile(path) || error("TLE file not found: $path")
    tles = read_tles_from_file(path)
    isempty(tles) && error("No valid TLE found in '$path'.")
    return tles
end

"""
    select_tles(tles, sel, path) -> Vector{TLE}

Filtra la lista de TLEs según el bloque `select` del grupo.

# Argumentos
- `tles::AbstractVector`: TLEs leídos del archivo.
- `sel`: bloque `select` (`norad_ids`, `names` o `range`), o `nothing` para usar todos.
- `path::AbstractString`: ruta del archivo, solo para los mensajes de error.

# Devuelve
- `Vector` con los TLEs elegidos, en el orden pedido (`norad_ids`, `names`) o del archivo (`range`).

# Excepciones
- `ErrorException`: si `select` tiene una clave no admitida, o no tiene exactamente una.
- `ErrorException`: si un NORAD o un nombre pedido no está en el archivo.
- `ErrorException`: si `range` no es `[primero, último]` válido dentro del archivo.

# Ejemplo
```julia
select_tles(tles, Dict("norad_ids" => [47699]), "ejemplo.tle")   # [TLE de AMAZONIA 1]
select_tles(tles, Dict("range" => [1, 50]), "catalogo.tle")       # los 50 primeros
```
"""
function select_tles(tles::AbstractVector, sel, path::AbstractString)
    sel === nothing && return collect(tles)

    modes = ("norad_ids", "names", "range")
    check_keys(sel, modes; who = "tle select")
    count(k -> haskey(sel, k), modes) == 1 ||
        error("'tle select': use exactly one of norad_ids, names or range.")

    if haskey(sel, "norad_ids")
        # Búsqueda por número NORAD, conservando el orden pedido
        return map(sel["norad_ids"]) do id
            i = findfirst(t -> t.satellite_number == Int(id), tles)
            i === nothing && error("NORAD id $id not found in '$path'.")
            tles[i]
        end
    elseif haskey(sel, "names")
        # Búsqueda por nombre exacto (sin espacios sobrantes)
        return map(sel["names"]) do name
            i = findfirst(t -> strip(t.name) == strip(String(name)), tles)
            i === nothing && error("TLE name '$name' not found in '$path'.")
            tles[i]
        end
    else
        r = sel["range"]
        length(r) == 2 || error("'tle select': 'range' must be [first, last].")
        a, b = Int(r[1]), Int(r[2])
        (1 <= a <= b <= length(tles)) ||
            error("'tle select': range [$a, $b] outside 1..$(length(tles)) (TLEs in '$path').")
        return tles[a:b]
    end
end

"""
    tle_to_elements(tle, ctx) -> Dict{String,Any}

Propaga un TLE con SGP4 hasta `ctx.epoch`, pasa el estado de TEME a J2000 y lo
convierte en elementos keplerianos osculadores.

Por qué existe: los elementos de un TLE son elementos MEDIOS del modelo SGP4, en
TEME y en la época del TLE; el simulador necesita elementos osculadores en J2000
y en la época común del escenario.

# Argumentos
- `tle::TLE`: TLE leído con `load_tles`.
- `ctx`: contexto del Builder (usa `ctx.epoch` y `ctx.mu`).

# Devuelve
- `Dict{String,Any}` con `SMA, ECC, INC, RAAN, AOP, TA` y `_label` = número NORAD
  (el Builder lo usa para el nombre y no lo escribe en el JSON).

# Excepciones
- `ErrorException`: si SGP4 falla (por ejemplo, satélite reentrado a esa fecha).
- `ErrorException`: si el estado resultante no es una órbita ligada.

# Ejemplo
```julia
ctx = (epoch = DateTime("2023-03-24T18:00:00"), mu = 398600.4415)
tle_to_elements(load_tles("config/tle/ejemplo.tle")[1], ctx)
# SMA ≈ 7135.25, ECC ≈ 0.00094, INC ≈ 98.470, RAAN ≈ 161.897, AOP + TA ≈ 328.92, _label = "47699"
```
"""
function tle_to_elements(tle::TLE, ctx)
    id = "$(tle.satellite_number) ($(strip(tle.name)))"          # identificador para mensajes

    # Tiempo desde la época del TLE hasta init_epoch, en minutos (puede ser negativo)
    Δt_min = Dates.value(ctx.epoch - tle.epoch) / 60_000
    abs(Δt_min) / 1440 > TLE_MAX_DAYS &&
        @warn "TLE $id: init_epoch is $(round(Δt_min / 1440, digits = 1)) days from the TLE epoch; SGP4 accuracy degrades."

    # 1) SGP4: posición y velocidad en TEME (km, km/s)
    r_teme, v_teme = try
        sgp4!(sgp4_init(tle), Δt_min)
    catch err
        error("TLE $id: SGP4 propagation failed: " * sprint(showerror, err))
    end

    # 2) TEME -> J2000 en init_epoch; la misma matriz vale para la velocidad
    #    (el giro de TEME por precesión-nutación es ~1e-12 rad/s, despreciable)
    D = r_eci_to_eci(TEME(), J2000(), datetime2julian(ctx.epoch))
    r = Vector{Float64}(D * r_teme)
    v = Vector{Float64}(D * v_teme)

    # 3) Elementos osculadores con el μ del simulador (es el que usará al propagar)
    el = cartesian_to_kepler(r, v, ctx.mu)
    el["_label"] = string(tle.satellite_number)
    return el
end

"""
    tle_pattern(group, ctx) -> Vector{Dict{String,Any}}

Genera un miembro por cada TLE seleccionado del archivo del grupo.

# Argumentos
- `group::AbstractDict`: grupo de la spec con `file` y, opcionalmente, `select`.
- `ctx`: contexto del Builder (usa `ctx.spec_dir`, `ctx.epoch`, `ctx.mu`).

# Devuelve
- `Vector` de diccionarios de elementos, uno por TLE seleccionado.

# Excepciones
- `ErrorException`: si falta `file`, el archivo no existe o la selección no es válida.
- `ErrorException`: cualquier error de `tle_to_elements`.

# Ejemplo
```julia
group = Dict("file" => "../tle/ejemplo.tle", "select" => Dict("norad_ids" => [47699]))
tle_pattern(group, ctx)        # 1 miembro; el Builder lo llamará <grupo>_47699
```
"""
function tle_pattern(group::AbstractDict, ctx)
    haskey(group, "file") || error("Pattern 'tle' needs a 'file' with the TLEs.")

    # Ruta relativa a la carpeta de la spec, igual que "template"
    path = String(group["file"])
    isabspath(path) || (path = normpath(joinpath(ctx.spec_dir, path)))

    tles = select_tles(load_tles(path), get(group, "select", nothing), path)
    return [tle_to_elements(t, ctx) for t in tles]
end

register_pattern!("tle", tle_pattern)
