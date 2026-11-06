# =============================================================================
# Visibility.jl
#
# Módulo Visibility: condiciones geométricas para que un objetivo se vea desde un
# observador en un instante. Las usa Photometry antes de calcular la magnitud.
#
# Funciones:
#   eclipseEarthECI : ¿está el objetivo en la sombra de la Tierra? (umbra / penumbra)
#   losECI          : ¿hay línea de visión entre dos puntos, sin que la Tierra la corte?
#   sunGlareConeECI : ¿está el objetivo dentro del cono de deslumbramiento del Sol?
#
# Modelo: Tierra esférica de radio R (Constants), Sol a distancia D_sun con radio R_sun,
# todo en ECI y en km. La dirección del Sol es la geocéntrica (la diferencia con la
# vista desde el satélite es < 0,003° en LEO-GEO).
#
# Dependencias: Constants (R, R_sun, D_sun).
# =============================================================================

module Visibility

export eclipseEarthECI, losECI, sunGlareConeECI

using ..Constants: R, R_sun, D_sun

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

    r_u(x) = R (1 - x / L_u),   L_u = D_sun R / (R_sun - R)
    r_p(x) = R (1 + x / L_p),   L_p = D_sun R / (R_sun + R)

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

    L_umbra = D_sun * R / (R_sun - R)                # vértice del cono de umbra
    L_penum = D_sun * R / (R_sun + R)                # vértice (delante de la Tierra) del de penumbra
    r_umbra = max(R * (1.0 - x / L_umbra), 0.0)
    r_penum = R * (1.0 + x / L_penum)

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
    depth = R - rn                                   # > 0 si está bajo la superficie
    depth <= TOL && return r
    depth <= MAX_PROJECT_DEPTH_KM && return (R / rn) * r
    return nothing
end

"""
    losECI(r1, r2) -> Bool

Indica si hay línea de visión entre dos puntos, es decir, si el segmento que los une
no atraviesa la Tierra (esfera de radio R). Con r(s) = r1 + s (r2 - r1), s ∈ [0, 1],
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
    sunGlareConeECI(u_sun, r_obs, r_obj, angMaxDeg) -> (Bool, Float64)

Indica si el objetivo, visto desde el observador, está a menos de `angMaxDeg` del Sol,
dentro del cono en el que la luz solar deslumbra la cámara.

Por qué existe: aunque el objetivo esté iluminado y en LOS, cerca del Sol la cámara no
puede observarlo.

# Argumentos
- `u_sun::SVector{3,Float64}`: dirección unitaria hacia el Sol en ECI.
- `r_obs::SVector{3,Float64}`, `r_obj::SVector{3,Float64}`: posiciones en ECI (km).
- `angMaxDeg::Real`: semiángulo del cono de deslumbramiento (°).

# Devuelve
- `(deslumbrado, angulo_deg)`: el ángulo Sol-objetivo visto desde el observador (°).

# Excepciones
- Ninguna (si `r_obj == r_obs` el ángulo sale `NaN` y no hay deslumbramiento).
"""
function sunGlareConeECI(u_sun::SVector{3,Float64}, r_obs::SVector{3,Float64},
                         r_obj::SVector{3,Float64}, angMaxDeg::Real)
    u_obj = normalize(r_obj - r_obs)                 # dirección observador -> objetivo
    ang_deg = acosd(clamp(dot(u_sun, u_obj), -1.0, 1.0))
    TOL_DEG = 1e-6                                   # tolerancia angular (°)
    return ang_deg <= angMaxDeg + TOL_DEG, ang_deg
end

end # module Visibility
