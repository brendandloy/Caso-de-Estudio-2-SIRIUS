# =============================================================================
# Visibility.jl
#
# Módulo Visibility: condiciones geométricas para que un objetivo se vea desde un
# observador en un instante (las usa Photometry antes de calcular la magnitud) y para
# que la cámara pueda observar (la usa CameraSim).
#
# Funciones:
#   shadow_function         : fracción ν del disco solar que ve un punto (1 iluminado, 0 umbra)
#   losECI                  : ¿hay línea de visión entre dos puntos, sin que la Tierra la corte?
#   sun_boresight_angle_deg : ángulo entre el eje óptico de la cámara y el Sol (exclusión solar)
#   earth_occulting_radius  : radio de la esfera que tapa la línea de visión (Tierra + atmósfera)
#   earth_angular_radius_deg: semiángulo del disco de la Tierra visto desde un punto
#   earth_ray_hit           : punto donde un rayo desde el observador choca con la Tierra
#   is_earth_occulted       : ¿tapa la Tierra una dirección (estrella) vista desde el observador?
#
# Modelo: Tierra esférica de radio R = R_EARTH, Sol a distancia D = D_SUN con radio R☉ = R_SUN,
# todo en ECI y en km. Para tapar la línea de visión se usa R_EARTH + EARTH_ATMOSPHERE_ALTITUDE_KM
# (la atmósfera), salvo que un extremo esté por debajo de esa altura (earth_occulting_radius). La dirección del Sol es la geocéntrica (la diferencia con la
# vista desde el satélite es < 0,003° en LEO-GEO).
#
# Dependencias: Constants (R_EARTH, R_SUN, D_SUN) y Parameters (EARTH_ATMOSPHERE_ALTITUDE_KM).
# =============================================================================

module Visibility

export shadow_function, losECI, sun_boresight_angle_deg
export earth_occulting_radius, earth_angular_radius_deg, earth_ray_hit, is_earth_occulted

using ..Constants: R_EARTH, R_SUN, D_SUN
using ..Parameters: EARTH_ATMOSPHERE_ALTITUDE_KM

using LinearAlgebra
using StaticArrays

# Tolerancia numérica común (km o adimensional, según el caso)
const TOL = 1e-9

# Profundidad máxima (km) bajo la superficie que se admite en un extremo de la LOS
# antes de proyectarlo sobre ella (p. ej. estaciones con altura elipsoidal negativa)
const MAX_PROJECT_DEPTH_KM = 100.0

"""
    shadow_function(r_obj, u_sun) -> Float64

Función de sombra ν de la Tierra: fracción del disco solar que ve un punto, de 1
(iluminado) a 0 (umbra), con valores intermedios en penumbra. Se modela como el solape de
dos discos vistos desde el punto (Montenbruck y Gill, 2000, sec. 3.4.2, ecs. 3.85–3.94):

    a = arcsin(R☉ / |r☉ - r|)                     radio aparente del Sol
    b = arcsin(R⊕ / |r|)                           radio aparente de la Tierra
    c = arccos(-r·(r☉ - r) / (|r| |r☉ - r|))       separación de sus centros

- c ≥ a + b: no hay ocultación, ν = 1.
- c ≤ b - a: ocultación total (umbra), ν = 0.
- c ≤ a - b: anular (la Tierra entera dentro del disco solar, solo a más de ~1,4·10⁶ km),
  ν = 1 - b²/a².
- Si no, parcial: con x = (c² + a² - b²)/(2c) e y = √(a² - x²), el área tapada es
  A = a² arccos(x/a) + b² arccos((c - x)/b) - c y, y ν = 1 - A/(π a²).

Se usa la Tierra sin atmósfera (radio R⊕), como Montenbruck y Gill, y el Sol en
r☉ = D☉ û☉.

Por qué existe: el flujo que refleja el objetivo es proporcional a la luz solar que le
llega, así que la fotometría lo multiplica por ν (la magnitud aumenta en -2,5 log₁₀ ν).

# Argumentos
- `r_obj::SVector{3,Float64}`: posición del punto en ECI (km).
- `u_sun::SVector{3,Float64}`: dirección unitaria Tierra -> Sol en ECI.

# Devuelve
- `Float64` ν ∈ [0, 1].

# Excepciones
- Ninguna.

# Ejemplo
```julia
shadow_function(SVector(-7000.0, 0, 0), SVector(1.0, 0, 0))   # 0.0: en el eje de la sombra (umbra)
```
"""
function shadow_function(r_obj::SVector{3,Float64}, u_sun::SVector{3,Float64})
    r_sun_rel = D_SUN * u_sun - r_obj                  # del punto al Sol
    d_sun, d_earth = norm(r_sun_rel), norm(r_obj)
    a = asin(R_SUN / d_sun)
    b = asin(min(R_EARTH / d_earth, 1.0))
    c = acos(clamp(-dot(r_obj, r_sun_rel) / (d_earth * d_sun), -1.0, 1.0))

    c >= a + b && return 1.0                           # sin ocultación
    c <= b - a && return 0.0                           # umbra
    c <= a - b && return 1.0 - (b / a)^2               # anular
    return 1.0 - occulted_area(a, b, c) / (π * a^2)
end

# Área del disco solar (radio aparente a) tapada por el de la Tierra (b) con sus centros a c,
# para |a - b| < c < a + b (Montenbruck y Gill, 2000, ecs. 3.92 y 3.93); ángulos en rad
function occulted_area(a::Real, b::Real, c::Real)
    x = (c^2 + a^2 - b^2) / (2c)
    y = sqrt(max(a^2 - x^2, 0.0))
    return a^2 * acos(clamp(x / a, -1.0, 1.0)) + b^2 * acos(clamp((c - x) / b, -1.0, 1.0)) - c * y
end

"""
    project_to_surface(r) -> Union{SVector{3,Float64}, Nothing}

Devuelve `r` si está sobre la superficie terrestre; si está como mucho
`MAX_PROJECT_DEPTH_KM` por debajo, su proyección radial sobre la esfera de radio R;
y `nothing` si está más hondo.

Por qué existe: un extremo justo bajo la superficie (una estación con altura
negativa) cortaría siempre la esfera y nunca tendría LOS.

# Argumentos
- `r::SVector{3,Float64}`: punto en ECI (km).

# Devuelve
- `SVector{3,Float64}` (el punto o su proyección), o `nothing`.

# Excepciones
- Ninguna.
"""
function project_to_surface(r::SVector{3,Float64})
    rn = norm(r)
    depth = R_EARTH - rn                                   # > 0 si está bajo la superficie
    depth <= TOL && return r
    depth <= MAX_PROJECT_DEPTH_KM && return (R_EARTH / rn) * r
    return nothing
end

"""
    losECI(r1, r2) -> Bool

Indica si hay línea de visión entre dos puntos, es decir, si el segmento que los une
no atraviesa la Tierra (esfera de radio R = `earth_occulting_radius(r1, r2)`: la Tierra
más la atmósfera si los dos puntos están por encima de ella). Con r(s) = r1 + s (r2 - r1), s ∈ [0, 1],
la condición |r(s)| = R da la ecuación de segundo grado a s² + b s + c = 0, con

    a = |r2 - r1|²,   b = 2 r1·(r2 - r1),   c = |r1|² - R².

Si no tiene raíces reales, o las raíces caen fuera de (0, 1), hay LOS. Una tangencia
dentro del segmento cuenta como bloqueada.

Por qué existe: un objetivo detrás de la Tierra no se ve aunque esté iluminado.

# Argumentos
- `r1::SVector{3,Float64}`, `r2::SVector{3,Float64}`: extremos en ECI (km). Un extremo
  hasta `MAX_PROJECT_DEPTH_KM` bajo la superficie se proyecta sobre ella.

# Devuelve
- `true` si hay LOS; `false` si la Tierra la corta, si un extremo está demasiado hondo
  o si los dos puntos coinciden.

# Excepciones
- Ninguna.

# Ejemplo
```julia
losECI(SVector(7000.0, 0, 0), SVector(-7000.0, 0, 0))   # false: la Tierra está en medio
```
"""
function losECI(r1::SVector{3,Float64}, r2::SVector{3,Float64})
    p1 = project_to_surface(r1)
    p2 = project_to_surface(r2)
    (p1 === nothing || p2 === nothing) && return false   # extremo demasiado bajo la superficie

    # Coeficientes de |p1 + s (p2 - p1)|² = R²
    R = earth_occulting_radius(p1, p2)
    d = p2 - p1
    a = dot(d, d)
    a <= TOL && return false                          # puntos coincidentes
    b = 2.0 * dot(p1, d)
    c = dot(p1, p1) - R^2
    disc = b^2 - 4.0 * a * c
    scale = b^2 + 4.0 * a * abs(c)                    # escala del discriminante, para las tolerancias

    disc < -scale * TOL && return true                # la recta no corta la esfera

    # Tangencia: bloquea solo si el punto de contacto está dentro del segmento
    if abs(disc) <= max(1.0, scale) * TOL
        s = -b / (2.0 * a)
        return !(TOL < s < 1.0 - TOL)
    end

    # Dos cortes: bloquea si el tramo [s1, s2] dentro de la esfera solapa con (0, 1)
    sq = sqrt(disc)
    s1, s2 = (-b - sq) / (2.0 * a), (-b + sq) / (2.0 * a)
    return !(max(s1, s2) > TOL && min(s1, s2) < 1.0 - TOL)
end

"""
    earth_occulting_radius(points...) -> Float64

Radio (km) de la esfera con la que la Tierra tapa la línea de visión desde los puntos dados
(ECI, km): `R_EARTH + EARTH_ATMOSPHERE_ALTITUDE_KM` si todos están por encima de la
atmósfera y `R_EARTH` si alguno está dentro.

Por qué existe: cerca del limbo la atmósfera absorbe y desvía la luz, así que tapa como si
la Tierra fuera más grande; pero un punto dentro de ella (una estación en tierra) tendría
entonces todo el cielo tapado, y para él se usa la superficie.
"""
function earth_occulting_radius(points::AbstractVector{<:Real}...)
    R_atm = R_EARTH + EARTH_ATMOSPHERE_ALTITUDE_KM
    return all(p -> norm(p) > R_atm, points) ? R_atm : R_EARTH
end

"""
    earth_angular_radius_deg(r_obs) -> Float64

Semiángulo (°) del disco de la Tierra visto desde `r_obs` (ECI, km), centrado en el nadir:
sin ρ⊕ = R / |r_obs|, con R = `earth_occulting_radius(r_obs)` (la tangente a la esfera).
Vale 90° si el punto está sobre la esfera o dentro.
"""
earth_angular_radius_deg(r_obs::AbstractVector{<:Real}) = asind(min(earth_occulting_radius(r_obs) / norm(r_obs), 1.0))

"""
    earth_ray_hit(r_obs, u) -> Union{SVector{3,Float64}, Nothing}

Primer punto donde el rayo que sale de `r_obs` en la dirección `u` choca con la esfera de
la Tierra (radio `earth_occulting_radius(r_obs)`), o `nothing` si no la toca. Con
r(s) = r_obs + s u y |u| = 1, |r(s)|² = R² da s² + 2b s + c = 0, con b = r_obs·u y
c = |r_obs|² - R²; hay choque delante del observador si c > 0, b < 0 y b² - c > 0, en
s = -b - √(b² - c). Una tangencia no cuenta como choque.

Por qué existe: decide qué estrellas tapa la Tierra (están en el infinito, así que no sirve
`losECI`) y en qué punto de la Tierra cae cada píxel, para saber si es de día.

# Argumentos
- `r_obs::AbstractVector{<:Real}`: posición del observador en ECI (km).
- `u::AbstractVector{<:Real}`: dirección unitaria del rayo en ECI.

# Devuelve
- El punto de choque en ECI (km), o `nothing` (también si el observador está bajo la superficie).

# Excepciones
- Ninguna.

# Ejemplo
```julia
earth_ray_hit(SVector(10000.0, 0, 0), SVector(-1.0, 0, 0))   # [R, 0, 0], con R = R_EARTH + 100
```
"""
function earth_ray_hit(r_obs::AbstractVector{<:Real}, u::AbstractVector{<:Real})
    R = earth_occulting_radius(r_obs)
    b = dot(r_obs, u)
    c = dot(r_obs, r_obs) - R^2
    c <= 0 && return nothing                  # observador bajo la superficie
    disc = b^2 - c
    (b >= 0 || disc <= 0) && return nothing   # la esfera queda detrás, o el rayo no la corta
    return SVector{3,Float64}(r_obs + (-b - sqrt(disc)) * u)
end

# ¿Tapa la Tierra (con la atmósfera) la dirección unitaria u (ECI) vista desde r_obs (ECI, km)?
is_earth_occulted(r_obs::AbstractVector{<:Real}, u::AbstractVector{<:Real}) = earth_ray_hit(r_obs, u) !== nothing

"""
    sun_boresight_angle_deg(C_ci, u_sun) -> Float64

Ángulo (°) entre el eje óptico (boresight) de la cámara y la dirección del Sol. El
boresight es z_c = [0, 0, 1] en ejes cámara; como C_ci pasa de ECI a cámara, en ECI es
C_ciᵀ [0, 0, 1]ᵀ, la tercera fila de C_ci, y el ángulo es θ☉ = arccos(b̂ · û☉).

Por qué existe: si θ☉ es menor que el semiángulo de exclusión del parasol de la cámara
(`camera_config.sun_exclusion_deg`), la luz del Sol entra en la óptica y satura toda la
imagen, así que la cámara no observa en ese instante.

# Argumentos
- `C_ci::AbstractMatrix`: DCM ECI -> cámara (3 × 3).
- `u_sun::SVector{3,Float64}`: dirección unitaria del Sol en ECI (geocéntrica; el paralaje
  desde la órbita es < 0,003° hasta GEO).

# Devuelve
- `Float64`: θ☉ en grados, en [0, 180].

# Excepciones
- Ninguna.

# Ejemplo
```julia
sun_boresight_angle_deg(I(3), SVector(0.0, 0.0, 1.0))   # 0.0: el Sol en el centro de la imagen
```
"""
function sun_boresight_angle_deg(C_ci::AbstractMatrix, u_sun::SVector{3,Float64})
    boresight = SVector(C_ci[3, 1], C_ci[3, 2], C_ci[3, 3])     # z_c en ECI: tercera fila de C_ci
    return acosd(clamp(dot(boresight, u_sun), -1.0, 1.0))
end

end # module Visibility
