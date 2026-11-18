# =============================================================================
# Visibility.jl
#
# Módulo Visibility: condiciones geométricas para que un objetivo se vea desde un
# observador en un instante (las usa Photometry antes de calcular la magnitud) y para
# que la cámara pueda observar (la usa CameraSim).
#
# Funciones:
#   eclipseEarthECI         : ¿está el objetivo en la sombra de la Tierra? (umbra / penumbra)
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

export eclipseEarthECI, losECI, sun_boresight_angle_deg
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
    eclipseEarthECI(r_sc, u_sun) -> (Bool, String)

Indica si un punto está en la sombra de la Tierra con el modelo de conos de umbra y
penumbra. En el eje de la sombra (opuesto al Sol), a una distancia x detrás de la
Tierra, los radios de los conos son

    r_u(x) = R (1 - x / L_u),   L_u = D R / (R☉ - R)
    r_p(x) = R (1 + x / L_p),   L_p = D R / (R☉ + R)

con L_u y L_p las distancias del centro de la Tierra al vértice de cada cono. El
punto está en umbra si su distancia al eje ρ ≤ r_u(x), y en penumbra si ρ ≤ r_p(x).

Por qué existe: un objetivo en sombra no refleja luz del Sol, así que no aparece en
la imagen (la fotometría le da magnitud infinita, también en penumbra).

# Argumentos
- `r_sc::SVector{3,Float64}`: posición del punto en ECI (km).
- `u_sun::SVector{3,Float64}`: dirección unitaria Tierra -> Sol en ECI.

# Devuelve
- `(en_sombra, region)`: `region` es `"umbra"`, `"penumbra"` o `"lit"`.

# Excepciones
- Ninguna.

# Ejemplo
```julia
eclipseEarthECI(SVector(-7000.0, 0, 0), SVector(1.0, 0, 0))   # (true, "umbra")
```
"""
function eclipseEarthECI(r_sc::SVector{3,Float64}, u_sun::SVector{3,Float64})
    u_axis = -u_sun                                  # eje de la sombra, alejándose del Sol
    x = dot(r_sc, u_axis)                            # distancia detrás de la Tierra a lo largo del eje
    rho = norm(r_sc - x * u_axis)                    # distancia al eje
    x <= 0.0 && return false, "lit"                  # en el lado del Sol

    L_umbra = D_SUN * R_EARTH / (R_SUN - R_EARTH)                # vértice del cono de umbra
    L_penum = D_SUN * R_EARTH / (R_SUN + R_EARTH)                # vértice (delante de la Tierra) del de penumbra
    r_umbra = max(R_EARTH * (1.0 - x / L_umbra), 0.0)
    r_penum = R_EARTH * (1.0 + x / L_penum)

    (x <= L_umbra + TOL && r_umbra > 0.0 && rho <= r_umbra + TOL) && return true, "umbra"
    rho <= r_penum + TOL && return true, "penumbra"
    return false, "lit"
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
