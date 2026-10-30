# =============================================================================
# Quaternions.jl
#
# Parametrizaciones de actitud del tema T6 (Attitude Kinematics) y sorteo de
# actitudes aleatorias. Se incluye desde Attitude.jl.
#
# Convención: todas las DCM de este archivo son del tipo C_bi del T6, es decir,
# pasan componentes del sistema de referencia (ECI o LVLH) al sistema cuerpo.
#
#   Matrices elementales (T6, diap. 19), rotación α del sistema alrededor de x, y, z:
#     C_x(α) = [1 0 0; 0 c s; 0 -s c]   C_y(α) = [c 0 -s; 0 1 0; s 0 c]   C_z(α) = [c s 0; -s c 0; 0 0 1]
#   Euler (T6, diap. 22): secuencia "abc", ángulos [θ1, θ2, θ3] en orden de aplicación
#     C = C_c(θ3) · C_b(θ2) · C_a(θ1)        p. ej. "321" -> C = C_x(φ) C_y(θ) C_z(ψ)
#   Cuaternión (T6, diap. 36-37): q = [cos(φ/2), e sin(φ/2)], escalar primero
#     C(q) = (q_η² - q_vᵀq_v) I + 2 q_v q_vᵀ - 2 q_η [q_v×]
#   Cinemática (T6, diap. 49), con ω en ejes cuerpo:
#     q̇ = ½ Ω(ω) q,   Ω(ω) = [0 -ωᵀ; ω -[ω×]]
#
# Funciones:
#   skew, quat_to_dcm, elementary_dcm, euler_to_dcm  : parametrizaciones -> DCM
#   dcm_to_quat                                     : DCM -> cuaternión (método de Shepperd)
#   omega_matrix, quat_step                         : propagación del cuaternión paso a paso
#   parse_orientation                               : lee "q" o "euler_deg" + "sequence" de un bloque
#   random_quaternion, random_unit_vector           : sorteo uniforme de actitudes y ejes
# =============================================================================

# Matriz antisimétrica [v×] tal que [v×] w = v × w
skew(v::AbstractVector) = @SMatrix [0.0 -v[3] v[2]; v[3] 0.0 -v[1]; -v[2] v[1] 0.0]

# Tolerancia para aceptar un cuaternión escrito a mano como unitario (luego se normaliza)
const QUAT_NORM_TOL = 1e-3

"""
    quat_to_dcm(q) -> SMatrix{3,3,Float64,9}

Convierte un cuaternión (escalar primero) en la DCM C_bi con la fórmula del T6:
`C(q) = (q_η² - q_vᵀq_v) I + 2 q_v q_vᵀ - 2 q_η [q_v×]`.

El cuaternión se normaliza antes de usarlo, así se admiten valores escritos a
mano con pocos decimales (p. ej. `[0.7071, 0, 0, 0.7071]`).

# Argumentos
- `q::AbstractVector`: `[q_η, q_1, q_2, q_3]`.

# Devuelve
- `SMatrix{3,3,Float64,9}`: C_bi (referencia -> cuerpo).

# Excepciones
- `ErrorException`: si `q` no tiene 4 componentes o su norma se aleja de 1 más de `QUAT_NORM_TOL`.

# Ejemplo
```julia
# 90° alrededor de z: el eje x del cuerpo apunta al eje y del sistema de referencia
C = quat_to_dcm([cosd(45), 0, 0, sind(45)])
C' * [1, 0, 0]          # ≈ [0, 1, 0]
```
"""
function quat_to_dcm(q::AbstractVector)
    length(q) == 4 || error("Quaternion must have 4 components [q0, q1, q2, q3] (got $(length(q))).")
    n = norm(q)
    abs(n - 1) <= QUAT_NORM_TOL || error("Quaternion is not unitary (norm = $(round(n, sigdigits = 6))).")

    qe = q[1] / n                                       # parte escalar q_η
    qv = SVector{3,Float64}(q[2], q[3], q[4]) / n       # parte vectorial q_v
    I3 = SMatrix{3,3,Float64,9}(I)

    return (qe^2 - dot(qv, qv)) * I3 + 2 * qv * qv' - 2 * qe * skew(qv)
end

"""
    dcm_to_quat(C) -> SVector{4,Float64}

Cuaternión (escalar primero, `q_η >= 0`) de una DCM C_bi. Generaliza la fórmula del
T6 (diap. 36), `q_η = ½√(1 + tr C)`, `q_v = [C23 - C32, C31 - C13, C12 - C21] / (4 q_η)`,
con el método de Shepperd: se despeja primero la componente de mayor módulo
(q_η, q1, q2 o q3) para no dividir por un número pequeño. La fórmula del T6 falla
cuando q_η ≈ 0, es decir, para giros de ≈ 180°.

# Argumentos
- `C::AbstractMatrix`: DCM 3 × 3 (C_bi, ortogonal y con determinante 1).

# Devuelve
- `SVector{4,Float64}` unitario con `q_η >= 0` (q y -q son la misma actitud).

# Excepciones
- Ninguna (no comprueba que `C` sea una rotación).

# Ejemplo
```julia
dcm_to_quat(elementary_dcm(3, π/2))      # ≈ [0.7071, 0, 0, 0.7071]
```
"""
function dcm_to_quat(C::AbstractMatrix)
    tr = C[1, 1] + C[2, 2] + C[3, 3]
    # 4·(componente)² de cada componente: se elige la mayor
    cand = (1 + tr, 1 + C[1, 1] - C[2, 2] - C[3, 3], 1 - C[1, 1] + C[2, 2] - C[3, 3], 1 - C[1, 1] - C[2, 2] + C[3, 3])
    k = argmax(cand)
    s = 0.5 * sqrt(cand[k])                               # componente de mayor módulo
    f = 0.25 / s

    q = if k == 1
        SVector(s, (C[2, 3] - C[3, 2]) * f, (C[3, 1] - C[1, 3]) * f, (C[1, 2] - C[2, 1]) * f)
    elseif k == 2
        SVector((C[2, 3] - C[3, 2]) * f, s, (C[1, 2] + C[2, 1]) * f, (C[1, 3] + C[3, 1]) * f)
    elseif k == 3
        SVector((C[3, 1] - C[1, 3]) * f, (C[1, 2] + C[2, 1]) * f, s, (C[2, 3] + C[3, 2]) * f)
    else
        SVector((C[1, 2] - C[2, 1]) * f, (C[1, 3] + C[3, 1]) * f, (C[2, 3] + C[3, 2]) * f, s)
    end

    q = q / norm(q)
    return q[1] < 0 ? -q : q                              # convenio q_η >= 0
end

"""
    omega_matrix(ω) -> SMatrix{4,4,Float64,16}

Matriz Ω(ω) de la ecuación cinemática del cuaternión (T6, diap. 49), para el
cuaternión escalar primero: `q̇ = ½ Ω(ω) q`, con `Ω = [0 -ωᵀ; ω -[ω×]]`.

# Argumentos
- `ω::AbstractVector`: velocidad angular del cuerpo respecto a la referencia, en ejes cuerpo (rad/s).

# Devuelve
- `SMatrix{4,4,Float64,16}`. Cumple `Ω² = -|ω|² I₄`.

# Excepciones
- Ninguna.
"""
function omega_matrix(ω::AbstractVector)
    wx, wy, wz = ω[1], ω[2], ω[3]
    return @SMatrix [0.0  -wx  -wy  -wz;
                     wx   0.0   wz  -wy;
                     wy  -wz   0.0   wx;
                     wz   wy  -wx   0.0]
end

"""
    quat_step(q, ω, Δt) -> SVector{4,Float64}

Avanza el cuaternión un paso de tiempo `Δt` con ω constante durante el paso:

    q(t + Δt) = exp(½ Ω(ω) Δt) · q(t)
              = [cos(|ω|Δt/2) I₄ + sin(|ω|Δt/2)/|ω| · Ω(ω)] · q(t)

Es la solución exacta de `q̇ = ½ Ω(ω) q` para ω constante: como `Ω² = -|ω|² I₄`, la
serie de la exponencial se separa en términos pares (coseno) e impares (seno).
Equivale a componer q con un giro de ángulo |ω|Δt alrededor de ω/|ω|.

Después se renormaliza: la operación conserva la norma en aritmética exacta, pero el
redondeo en coma flotante (~1e-16 por paso) se acumularía paso a paso.

Si en el futuro ω varía (dinámica con inercia), la misma función sirve con el ω de
cada paso; entonces es una aproximación de primer orden en `Δt`.

# Argumentos
- `q::AbstractVector`: cuaternión actual (C_bi, escalar primero).
- `ω::AbstractVector`: velocidad angular en ejes cuerpo durante el paso (rad/s).
- `Δt::Real`: paso de tiempo (s).

# Devuelve
- `SVector{4,Float64}` unitario.

# Excepciones
- Ninguna.

# Ejemplo
```julia
q1 = quat_step(SVector(1.0, 0, 0, 0), SVector(0, 0, deg2rad(90)), 1.0)    # 90° alrededor de z
quat_to_dcm(q1)' * [1, 0, 0]                                              # ≈ [0, 1, 0]
```
"""
function quat_step(q::AbstractVector, ω::AbstractVector, Δt::Real)
    q4 = SVector{4,Float64}(q)
    wn = norm(ω)
    iszero(wn) && return q4                               # sin giro: no cambia

    a = 0.5 * wn * Δt                                     # mitad del ángulo girado en el paso
    E = cos(a) * SMatrix{4,4,Float64,16}(I) + (sin(a) / wn) * omega_matrix(ω)
    qn = E * q4
    return qn / norm(qn)                                  # renormalización
end

"""
    elementary_dcm(axis, α) -> SMatrix{3,3,Float64,9}

Matriz elemental de cambio de base del T6 (diap. 19): rotación de ángulo `α`
del sistema alrededor del eje 1 (x), 2 (y) o 3 (z).

# Argumentos
- `axis::Int`: 1, 2 o 3.
- `α::Real`: ángulo en radianes.

# Devuelve
- `SMatrix{3,3,Float64,9}`: C_x(α), C_y(α) o C_z(α).

# Excepciones
- `ErrorException`: si `axis` no es 1, 2 o 3.

# Ejemplo
```julia
elementary_dcm(3, π/2) * [0, 1, 0]       # ≈ [1, 0, 0]: el eje y de partida, visto desde el sistema girado
```
"""
function elementary_dcm(axis::Int, α::Real)
    c, s = cos(α), sin(α)
    axis == 1 && return @SMatrix [1.0 0.0 0.0; 0.0 c s; 0.0 -s c]
    axis == 2 && return @SMatrix [c 0.0 -s; 0.0 1.0 0.0; s 0.0 c]
    axis == 3 && return @SMatrix [c s 0.0; -s c 0.0; 0.0 0.0 1.0]
    error("Elementary rotation axis must be 1, 2 or 3 (got $axis).")
end

"""
    euler_to_dcm(angles, sequence) -> SMatrix{3,3,Float64,9}

DCM C_bi a partir de tres ángulos de Euler, siguiendo el T6 (diap. 22): para la
secuencia `"abc"` y ángulos `[θ1, θ2, θ3]` (en el orden en que se aplican),
`C = C_c(θ3) · C_b(θ2) · C_a(θ1)`.

# Argumentos
- `angles::AbstractVector`: `[θ1, θ2, θ3]` en radianes, en orden de aplicación.
- `sequence::AbstractString`: tres dígitos de {1, 2, 3}, sin dos seguidos iguales ("321", "313"...).

# Devuelve
- `SMatrix{3,3,Float64,9}`: C_bi.

# Excepciones
- `ErrorException`: si la secuencia no es válida o no hay 3 ángulos.

# Ejemplo
```julia
# 3-2-1 (yaw ψ, pitch θ, roll φ), como en la diap. 26 del T6
C = euler_to_dcm(deg2rad.([30, 10, 5]), "321")      # = C_x(5°) C_y(10°) C_z(30°)
```
"""
function euler_to_dcm(angles::AbstractVector, sequence::AbstractString)
    length(angles) == 3 || error("Euler angles: 3 angles needed (got $(length(angles))).")
    ok = length(sequence) == 3 && all(c -> c in ('1', '2', '3'), sequence) &&
         sequence[1] != sequence[2] && sequence[2] != sequence[3]
    ok || error("Euler sequence must be 3 digits of 1/2/3 with no two consecutive equal, e.g. \"321\" or \"313\" (got \"$sequence\").")

    ax = [parse(Int, c) for c in sequence]          # ejes en orden de aplicación
    # La primera rotación queda a la derecha: C = C_c(θ3) C_b(θ2) C_a(θ1)
    return elementary_dcm(ax[3], angles[3]) * elementary_dcm(ax[2], angles[2]) * elementary_dcm(ax[1], angles[1])
end

"""
    parse_orientation(b, who; required) -> Union{SMatrix{3,3,Float64,9}, Nothing}

Lee una orientación de un bloque, escrita como cuaternión (`q`) o como ángulos
de Euler (`euler_deg` + `sequence`). Es la orientación del cuerpo respecto al
sistema de referencia del modo (ECI o LVLH), en las convenciones del T6.

# Argumentos
- `b::AbstractDict`: bloque leído del JSON.
- `who::AbstractString`: texto que identifica el bloque en el mensaje de error.
- `required::Bool`: si es `true`, falta de orientación es un error; si no, devuelve `nothing`.

# Devuelve
- `SMatrix{3,3,Float64,9}`: C_bi, o `nothing` si no hay orientación y no es obligatoria.

# Excepciones
- `ErrorException`: si aparecen `q` y `euler_deg` a la vez, si falta `sequence`
  (o aparece sin `euler_deg`), o si los valores no son válidos.

# Ejemplo
```julia
parse_orientation(Dict("euler_deg" => [90, 0, 0], "sequence" => "321"), "deb"; required = true)
parse_orientation(Dict("q" => [1, 0, 0, 0]), "deb"; required = true)      # identidad
```
"""
function parse_orientation(b::AbstractDict, who::AbstractString; required::Bool)
    has_q     = haskey(b, "q")
    has_euler = haskey(b, "euler_deg")

    has_q && has_euler && error("$who: give the orientation either as 'q' or as 'euler_deg', not both.")
    !has_euler && haskey(b, "sequence") && error("$who: 'sequence' without 'euler_deg'.")

    # Los errores de conversión se repiten con el nombre del objetivo delante
    with_who(f) = try
        f()
    catch err
        err isa ErrorException ? error("$who: " * err.msg) : rethrow()
    end

    if has_q
        q = getvec(b, "q", 4, who)
        return with_who(() -> quat_to_dcm(q))
    elseif has_euler
        haskey(b, "sequence") || error("$who: 'euler_deg' needs a 'sequence' (e.g. \"321\").")
        seq = b["sequence"]
        seq isa AbstractString || error("$who: 'sequence' must be a text such as \"321\" (got $(repr(seq))).")
        angles = deg2rad.(getvec(b, "euler_deg", 3, who))
        return with_who(() -> euler_to_dcm(angles, seq))
    end

    required && error("$who: missing orientation: 'q' = [q0, q1, q2, q3] or 'euler_deg' = [θ1, θ2, θ3] + 'sequence'.")
    return nothing
end

"""
    random_quaternion(rng) -> SVector{4,Float64}

Sortea un cuaternión unitario con distribución uniforme sobre todas las
orientaciones posibles.

Método: cuatro normales N(0,1) independientes forman un vector de R⁴ cuya
dirección es uniforme sobre la esfera unidad S³ (la densidad gaussiana solo
depende de la norma). Como los cuaterniones unitarios cubren las rotaciones
(q y -q son la misma), esto da una actitud uniforme. Sortear tres ángulos de
Euler uniformes NO lo hace: acumula orientaciones cerca de los polos.

# Argumentos
- `rng::AbstractRNG`: generador aleatorio (con semilla, para poder repetir).

# Devuelve
- `SVector{4,Float64}` unitario, escalar primero, con `q_η >= 0`.

# Excepciones
- Ninguna.

# Ejemplo
```julia
q = random_quaternion(MersenneTwister(42))
```
"""
function random_quaternion(rng::AbstractRNG)
    q = SVector{4,Float64}(randn(rng), randn(rng), randn(rng), randn(rng))
    q = q / norm(q)
    return q[1] < 0 ? -q : q              # q y -q son la misma actitud; se elige q_η >= 0
end

"""
    random_unit_vector(rng) -> SVector{3,Float64}

Sortea una dirección uniforme sobre la esfera unidad (mismo método que
`random_quaternion`, con tres normales).

# Argumentos
- `rng::AbstractRNG`: generador aleatorio.

# Devuelve
- `SVector{3,Float64}` unitario.

# Excepciones
- Ninguna.
"""
function random_unit_vector(rng::AbstractRNG)
    v = SVector{3,Float64}(randn(rng), randn(rng), randn(rng))
    return v / norm(v)
end
