module Environment

export load_config, parse_arguments, nav_err!, load_kernels, gs_positions, sun_dir_eci, export_data_csvs!

using ..Types
using ..Transformations
using ..Constants
using ..Attitude: AttitudeSpec, parse_member_attitude    # lectura del bloque "attitude" de cada miembro
using ..Attitude: camera_mount_block, parse_camera_mount # montaje de la cámara de cada observador (C_cb)
using ..Geometry: ShapeSpec, parse_member_shape          # lectura del bloque "shape" de cada objetivo

using StaticArrays
using SPICE
using Dates
using LinearAlgebra
using JSON3
using CSV
using DataFrames

#------------------------------------FUNCTIONS----------------------------------------------

const _KERNELS_LOADED = Ref(false)

function load_kernels(kernel_dir::AbstractString)
    if _KERNELS_LOADED[]
        return
    end

    furnsh(joinpath(kernel_dir, "de440s.bsp"))
    furnsh(joinpath(kernel_dir, "naif0012.tls.txt"))

    _KERNELS_LOADED[] = true
    return nothing
end


"""
    parse_arguments(args=ARGS) -> NamedTuple

Lee los argumentos de `scripts/main.jl`: la spec es obligatoria y la carpeta de
resultados, opcional.

# Argumentos
- `args`: `[spec.json]` o `[spec.json, output_path]`.

# Devuelve
- `(spec_path, output_path)`: rutas absolutas; sin `output_path`, `nothing` (main.jl usa entonces
  `output/<nombre de la spec>/`).

# Excepciones
- `ErrorException`: si no hay argumentos o hay más de dos (el mensaje indica el uso).

# Ejemplo
```julia
parse_arguments(["config/specs/ejemplo_leo.json"])     # output_path = nothing -> output/ejemplo_leo/
```
"""
function parse_arguments(args = ARGS)
    # La spec es obligatoria: así siempre queda claro qué escenario se simula
    1 <= length(args) <= 2 || error(
        "Usage: julia --project=. scripts/main.jl <spec.json> [<output_path>]"
    )

    return (
        spec_path   = abspath(args[1]),
        output_path = length(args) == 2 ? abspath(args[2]) : nothing
    )
end


"""
    load_config(path) -> NamedTuple

Lee el JSON completo de un escenario (el que escribe el generador) y lo convierte
en lo que necesita la simulación: estados iniciales en ECI, actitud, forma y
montaje de la cámara de cada miembro, estaciones y configuración de cámara y
fotometría. Imprime un resumen.

Por qué existe: es el único punto en el que el simulador lee el JSON; el resto de
módulos trabaja con los tipos ya validados.

# Argumentos
- `path::AbstractString`: JSON completo (normalmente en `config/generated/`).

# Devuelve
- `NamedTuple` con los tiempos (`Duration_hours`, `Data_time_step_seconds`,
  `init_epoch`), los flags, `camera_config`, `photometry_config`, `objects`
  (estado inicial de cada miembro), `obj_keys`, `obj_attitude`, `obj_shape`,
  `obs_keys`, `obs_attitude`, `obs_mount`, `obs_pos_sigma`, `obs_vel_sigma` y
  `stations_llh`.

# Excepciones
- Las de JSON3 (JSON mal formado) y las de validación de cada bloque
  (`parse_member_attitude`, `parse_member_shape`, `parse_camera_mount`).
"""
function load_config(path::AbstractString)

    # Se lee el texto en memoria: JSON3.read(path) mapea el archivo (mmap) y, mientras la sesión
    # de Julia siga abierta, Windows no deja sobrescribirlo (p. ej. al regenerar la spec)
    cfg = JSON3.read(read(path, String))

    objects = Dict{Symbol, State}()                 # estado inicial de objetivos y observadores (ECI)
    objectives = parse_objectives!(objects, cfg.objectives)
    observers = parse_observers!(objects, cfg.observers)
    stations_llh = Dict{Symbol, GroundStationLLH}(
        Symbol(name) => GroundStationLLH(s.lat, s.lon, s.h) for (name, s) in cfg.stations_llh)

    env = (
        Duration_hours         = cfg.Duration_hours,
        Data_time_step_seconds = cfg.Data_time_step_seconds,
        init_epoch             = DateTime(cfg.init_epoch),
        navigation_error_flag  = cfg.navigation_error_flag,
        graphic_display_flag   = cfg.graphic_display_flag,
        camera_config          = camera_settings(cfg),
        photometry_config      = cfg.photometry_config,
        objects                = objects,
        obs_attitude           = observers.attitude,
        obs_mount              = observers.mount,
        obs_pos_sigma          = observers.pos_sigma,
        obs_vel_sigma          = observers.vel_sigma,
        obs_keys               = observers.keys,
        obj_keys               = objectives.keys,
        obj_attitude           = objectives.attitude,
        obj_shape              = objectives.shape,
        stations_llh           = stations_llh
    )
    print_config_summary(env)
    return env
end

# camera_config del JSON más los flags que usa la cámara (están en la raíz del JSON)
camera_settings(cfg) = merge(NamedTuple(cfg.camera_config), (
    star_csv_flag                    = Bool(cfg.star_csv_flag),
    navigation_error_flag            = Bool(cfg.navigation_error_flag),
    ignore_fov_for_measurements_flag = Bool(cfg.ignore_fov_for_measurements_flag),
    image_generation_flag            = Bool(cfg.image_generation_flag)))

"""
    parse_objectives!(objects, entries) -> NamedTuple

Lee la sección `objectives` del JSON: añade a `objects` el estado inicial de cada
objetivo (clave `obj_<nombre>`) y valida sus bloques `attitude` y `shape`.

Por qué existe: separa la lectura de cada sección de `load_config`.

# Argumentos
- `objects::Dict{Symbol,State}`: estados iniciales; se modifica.
- `entries`: sección `objectives` del JSON (JSON3).

# Devuelve
- `(keys, attitude, shape)`: claves en el orden del JSON y `Dict` clave -> `AttitudeSpec` / `ShapeSpec`.

# Excepciones
- Las de `parse_member_attitude` y `parse_member_shape` (bloques ausentes o no válidos).
"""
function parse_objectives!(objects::Dict{Symbol,State}, entries)
    member_keys = Symbol[]                          # claves en el orden del JSON
    attitude = Dict{Symbol, AttitudeSpec}()         # actitud de cada objetivo (modo + parámetros)
    shape = Dict{Symbol, ShapeSpec}()               # forma de cada objetivo (tipo + medidas + semilla)

    for (name, dict) in entries
        key = Symbol("obj_$name")
        objects[key] = kepler_to_cartesian(dict)
        attitude[key] = parse_member_attitude(dict, String(name))
        shape[key] = parse_member_shape(dict, String(name))   # si aún hay algo por sortear, se sortea aquí
        push!(member_keys, key)
    end
    return (keys = member_keys, attitude = attitude, shape = shape)
end

"""
    parse_observers!(objects, entries) -> NamedTuple

Lee la sección `observers` del JSON: añade a `objects` el estado inicial de cada
observador (clave `obs_<nombre>`) y valida su actitud, el montaje de su cámara y
sus errores de navegación.

Por qué existe: separa la lectura de cada sección de `load_config`.

# Argumentos
- `objects::Dict{Symbol,State}`: estados iniciales; se modifica.
- `entries`: sección `observers` del JSON (JSON3).

# Devuelve
- `(keys, attitude, mount, pos_sigma, vel_sigma)`: claves en el orden del JSON y
  `Dict` clave -> `AttitudeSpec`, C_cb (cuerpo -> cámara) y σ de posición (km) y velocidad (km/s).

# Excepciones
- Las de `parse_member_attitude` y `parse_camera_mount`.
"""
function parse_observers!(objects::Dict{Symbol,State}, entries)
    member_keys = Symbol[]                          # claves en el orden del JSON
    attitude = Dict{Symbol, AttitudeSpec}()             # actitud del satélite de cada observador
    mount = Dict{Symbol, SMatrix{3,3,Float64,9}}()      # montaje de su cámara, C_cb (cuerpo -> cámara)
    pos_sigma = Dict{Symbol, SVector{3,Float64}}()
    vel_sigma = Dict{Symbol, SVector{3,Float64}}()

    for (name, dict) in entries
        key = Symbol("obs_$name")
        objects[key] = kepler_to_cartesian(dict)
        attitude[key] = parse_member_attitude(dict, String(name); role = "Observer")
        # "camera_mount" opcional: por defecto C_cb = I, con aviso
        mount[key] = parse_camera_mount(camera_mount_block(dict, String(name)), String(name))
        pos_sigma[key] = SVector{3,Float64}(dict.position_sigma_km)
        vel_sigma[key] = SVector{3,Float64}(dict.velocity_sigma_km_s)
        push!(member_keys, key)
    end
    return (keys = member_keys, attitude = attitude, mount = mount, pos_sigma = pos_sigma, vel_sigma = vel_sigma)
end

# Resumen de la configuración leída, en la terminal
print_config_summary(env) = println("""
    Loaded config:
    - Duration (h): $(env.Duration_hours)h
    - Step (s): $(env.Data_time_step_seconds)s
    - Epoch: $(env.init_epoch)
    - Navigation error: $(env.navigation_error_flag)
    - Graphics: $(env.graphic_display_flag)
    - Objects: $(keys(env.objects))
    - Stations: $(keys(env.stations_llh))
    """)


function nav_err!(objects::Dict{Symbol,State}, obs_keys::Vector{Symbol},
                  obs_pos_sigma::Dict{Symbol,SVector{3,Float64}}, obs_vel_sigma::Dict{Symbol,SVector{3,Float64}})

    for obs_key in obs_keys
        st = objects[obs_key]

        σr = obs_pos_sigma[obs_key]
        σv = obs_vel_sigma[obs_key]

        δr = σr .* SVector{3,Float64}(randn(3))
        δv = σv .* SVector{3,Float64}(randn(3))

        nav_name = Symbol("nav_", String(obs_key))

        objects[nav_name] = State(
            st.r + δr,
            st.v + δv
        )
    end

    return nothing
end


function gs_positions(stations::Dict{Symbol,GroundStationLLH}, start_time, tvec::AbstractVector{<:Real})
    out = Dict{Symbol, TrajectoryData}()
    t_out = collect(Float64.(tvec))

    # Compute θ₀ for the start time
    jd = datetime2julian(start_time)
    θ₀ = theta0_from_utc(jd)

    for (key, st) in stations
        # 1) Convert LLH → ECEF (one time per station)
        lat = deg2rad(st.lat_deg)
        lon = deg2rad(st.lon_deg)
        hkm = st.h_m / 1000.0
        r_ecef = llh_to_ecef(lat, lon, hkm)

        positions = rotate_z.(Ref(r_ecef), θ_g.(θ₀, t_out))

        # Velocities = zero (unless you want real ECI velocity later)
        vzero = SVector{3,Float64}(0.0, 0.0, 0.0)
        velocities = fill(vzero, length(t_out))

        # 3×N matrix (each column is one timestamp)
        out[key] = TrajectoryData(positions, velocities, t_out)
    end

    return out
end


function sun_dir_eci(t0::Float64, t_sim::Float64)
    et = t0 + t_sim
    state, _ = spkezr("SUN", et, "J2000", "NONE", "EARTH")
    v = state[1:3] / norm(state[1:3])
    return SVector{3,Float64}(v)
end


# Componente i (1, 2 o 3) de una lista de vectores 3D, como columna
component(vs, i) = [v[i] for v in vs]

"""
    build_rso_state_df(data, obj_keys) -> DataFrame

Tabla con el estado verdadero (ECI, km y km/s) de todos los objetivos, objetivo
tras objetivo: columnas `t`, `obj_name`, `r_x`…`r_z`, `v_x`…`v_z`.

Por qué existe: es `rso_states.csv`, la referencia de los objetivos para el Kalman.

# Argumentos
- `data::Dict{Symbol,TrajectoryData}`: trayectorias propagadas.
- `obj_keys::Vector{Symbol}`: objetivos, en el orden en que se escriben.

# Devuelve
- `DataFrame`.

# Excepciones
- `KeyError`: si falta la trayectoria de algún objetivo.
"""
function build_rso_state_df(data, obj_keys)
    trajs = [data[k] for k in obj_keys]
    r = vcat(SVector{3,Float64}[], (tr.r for tr in trajs)...)     # posiciones de todos, en orden
    v = vcat(SVector{3,Float64}[], (tr.v for tr in trajs)...)

    return DataFrame(
        t        = vcat(Float64[], (tr.t for tr in trajs)...),
        obj_name = vcat(String[], (fill(String(k), length(data[k].t)) for k in obj_keys)...),
        r_x = component(r, 1), r_y = component(r, 2), r_z = component(r, 3),
        v_x = component(v, 1), v_y = component(v, 2), v_z = component(v, 3)
    )
end

"""
    build_state_df(data, obs_name; has_nav_error) -> DataFrame

Tabla con el estado verdadero del observador (`*_true`) y, si `has_nav_error`, el estado
con error de navegación (`*_meas`, trayectoria `nav_<observador>`). ECI, km y km/s.

Por qué existe: es `state_measurements.csv` de cada observador.

# Argumentos
- `data::Dict{Symbol,TrajectoryData}`: trayectorias propagadas.
- `obs_name::Symbol`: observador.
- `has_nav_error::Bool`: añadir las columnas con error de navegación.

# Devuelve
- `DataFrame`.

# Excepciones
- `KeyError`: si `has_nav_error` es `true` y no existe la trayectoria `nav_<observador>`.
"""
function build_state_df(data, obs_name; has_nav_error::Bool)
    traj = data[obs_name]
    df = DataFrame(
        t = traj.t,
        r_x_true = component(traj.r, 1), r_y_true = component(traj.r, 2), r_z_true = component(traj.r, 3),
        v_x_true = component(traj.v, 1), v_y_true = component(traj.v, 2), v_z_true = component(traj.v, 3)
    )

    if has_nav_error
        nav = data[Symbol("nav_", String(obs_name))]           # misma órbita con error de navegación
        df.r_x_meas = component(nav.r, 1); df.r_y_meas = component(nav.r, 2); df.r_z_meas = component(nav.r, 3)
        df.v_x_meas = component(nav.v, 1); df.v_y_meas = component(nav.v, 2); df.v_z_meas = component(nav.v, 3)
    end
    return df
end

"""
    build_los_df(los_meas, obs_name; has_nav_error) -> DataFrame

Tabla con las medidas LOS de un observador, objetivo tras objetivo: dirección
verdadera (ECI, unitaria), magnitud y, si `has_nav_error`, la dirección medida con error.

Por qué existe: es `los_measurements.csv` de cada observador, la entrada del Kalman.

# Argumentos
- `los_meas`: medidas de todos los observadores (`camera_sim!`).
- `obs_name::Symbol`: observador.
- `has_nav_error::Bool`: añadir las columnas `dir_*_meas`.

# Devuelve
- `DataFrame` (vacío, con las columnas, si el observador no tiene medidas).

# Excepciones
- Ninguna.
"""
function build_los_df(los_meas, obs_name; has_nav_error::Bool)
    meas = haskey(los_meas, obs_name) ? los_meas[obs_name] : Dict{Symbol, Vector{LOSMeasurement}}()
    names = String[String(obj) for (obj, list) in meas for _ in list]
    ms = LOSMeasurement[m for (_, list) in meas for m in list]
    dir_true = [m.dir_true for m in ms]

    df = DataFrame(
        t = Float64[m.t for m in ms],
        obj_name = names,
        dir_x_true = component(dir_true, 1), dir_y_true = component(dir_true, 2), dir_z_true = component(dir_true, 3),
        Vmag = Float64[m.Vmag for m in ms]
    )

    if has_nav_error
        dir_meas = [m.dir_meas for m in ms]
        df.dir_x_meas = component(dir_meas, 1); df.dir_y_meas = component(dir_meas, 2); df.dir_z_meas = component(dir_meas, 3)
    end
    return df
end


function export_data_csvs!(data, los_meas, env; output_root=joinpath(@__DIR__, "..", "output"))

    has_nav_error = env.navigation_error_flag

    # Common RSO truth file
    rso_df = build_rso_state_df(data, env.obj_keys)
    CSV.write(joinpath(output_root, "rso_states.csv"), rso_df; delim=';')

    for obs_name in env.obs_keys
        outdir = joinpath(output_root, String(obs_name))
        mkpath(outdir)

        state_df = build_state_df(data, obs_name; has_nav_error=has_nav_error)
        CSV.write(joinpath(outdir, "state_measurements.csv"), state_df; delim=';')

        los_df = build_los_df(los_meas, obs_name; has_nav_error=has_nav_error)
        CSV.write(joinpath(outdir, "los_measurements.csv"), los_df; delim=';')
    end

    return nothing
end


end # module
