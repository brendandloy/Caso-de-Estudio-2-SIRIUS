# =============================================================================
# random.jl
#
# Modo "random": actitud inicial y giro aleatorios, típico de basura espacial.
# Sortea las condiciones iniciales y luego se propaga igual que "spin".
#
#   "attitude": {"type": "random", "seed": 42}
#
# - El rango de |ω| no lo elige el usuario: es el parámetro maestro
#   RANDOM_W_RANGE_DEG_S (src/Parameters.jl).
# - "seed" (opcional): semilla del sorteo. Sin ella se sortea una nueva y se imprime
#   para poder repetir la simulación. El generador de JSON pone una distinta por miembro.
#
# Qué se sortea (con un MersenneTwister(seed) propio de este objetivo):
#   1. Orientación inicial: cuaternión uniforme sobre todas las actitudes (random_quaternion).
#   2. Eje de giro en ejes cuerpo: dirección uniforme sobre la esfera (random_unit_vector).
#   3. |ω|: uniforme en RANDOM_W_RANGE_DEG_S.
# El sorteo (sample_random_state) está separado de la propagación (spin_history) para
# que un futuro modo con ecuaciones de Euler pueda reutilizarlo.
#
# Funciones: parse_random, sample_random_state, compute_attitude(::RandomSpec, ...)
# =============================================================================

"""
    RandomSpec(w_min_rad_s, w_max_rad_s, seed)

Parámetros validados del modo `random`.

# Campos
- `w_min_rad_s::Float64`, `w_max_rad_s::Float64`: rango del módulo de ω (rad/s), tomado de
  `RANDOM_W_RANGE_DEG_S` al leer el bloque.
- `seed::Int`: semilla del sorteo de este objetivo.
"""
struct RandomSpec <: AttitudeSpec
    w_min_rad_s::Float64
    w_max_rad_s::Float64
    seed::Int
end

# Claves admitidas en el bloque "random"
const RANDOM_KEYS = ("type", "seed")

"""
    parse_random(b, name) -> RandomSpec

Lee y valida un bloque `attitude` de tipo `random`. El rango de |ω| sale del
parámetro maestro `RANDOM_W_RANGE_DEG_S`. Si no hay `seed`, sortea una y la imprime
junto al nombre del objetivo.

# Argumentos
- `b::Dict{String,Any}`: bloque `attitude`.
- `name::AbstractString`: nombre del objetivo, para los mensajes.

# Devuelve
- `RandomSpec`.

# Excepciones
- `ErrorException`: si hay claves no admitidas (p. ej. `w_range_deg_s`, que ya no se
  da en el JSON) o `seed` no es un entero >= 0.

# Ejemplo
```julia
parse_random(Dict{String,Any}("type" => "random", "seed" => 42), "deb")
```
"""
function parse_random(b::Dict{String,Any}, name::AbstractString)
    who = "Objective '$name' (attitude 'random')"
    check_keys(b, RANDOM_KEYS, who)

    wmin, wmax = RANDOM_W_RANGE_DEG_S                     # parámetro maestro (°/s)

    if haskey(b, "seed")
        s = b["seed"]
        (s isa Real && !(s isa Bool) && isinteger(s) && s >= 0) ||
            error("$who: 'seed' must be an integer >= 0 (got $(repr(s))).")
        seed = Int(s)
    else
        # Sin semilla: se sortea una y se imprime para poder repetir este caso
        seed = Int(rand(RandomDevice(), UInt32))
        println("Objective '$name': random attitude without 'seed'; using seed = $seed ",
                "(add \"seed\": $seed to its 'attitude' block to repeat it).")
    end

    return RandomSpec(deg2rad(wmin), deg2rad(wmax), seed)
end

"""
    sample_random_state(spec) -> NamedTuple

Sortea las condiciones iniciales del modo `random` con la semilla del objetivo.

# Argumentos
- `spec::RandomSpec`: parámetros del modo.

# Devuelve
- `(q0, axis_b, w_rad_s)`: cuaternión inicial sorteado (C_bi, escalar primero), eje de
  giro unitario en ejes cuerpo y módulo de ω (rad/s).

# Excepciones
- Ninguna.

# Ejemplo
```julia
st = sample_random_state(RandomSpec(0.0, deg2rad(5), 42))
rad2deg(st.w_rad_s)        # entre 0 y 5 °/s; siempre el mismo con seed = 42
```
"""
function sample_random_state(spec::RandomSpec)
    rng = MersenneTwister(spec.seed)

    q0     = random_quaternion(rng)                      # orientación uniforme (C_bi)
    axis_b = random_unit_vector(rng)                     # eje de giro uniforme en ejes cuerpo
    w      = spec.w_min_rad_s + (spec.w_max_rad_s - spec.w_min_rad_s) * rand(rng)   # |ω| uniforme

    return (q0 = q0, axis_b = axis_b, w_rad_s = w)
end

"""
    compute_attitude(spec::RandomSpec, traj, sun_dir) -> AttitudeData

Actitud del modo `random`: sortea las condiciones iniciales y propaga el
cuaternión con ω constante (igual que `spin`, con `spin_history`).

# Argumentos
- `spec::RandomSpec`: parámetros del modo.
- `traj::TrajectoryData`: trayectoria del objetivo (solo se usan los tiempos).
- `sun_dir::AbstractVector`: dirección del Sol (no se usa).

# Devuelve
- `AttitudeData`.

# Excepciones
- Ninguna.
"""
function compute_attitude(spec::RandomSpec, traj::TrajectoryData, sun_dir::AbstractVector)
    st = sample_random_state(spec)
    return AttitudeData(spin_history(st.q0, st.axis_b, st.w_rad_s, traj.t))
end

register_attitude!("random", parse_random)
