# =============================================================================
# Photometry.jl
#
# Módulo Photometry: magnitud aparente de cada objetivo visto desde cada observador
# en cada instante, con el modelo de reflexión de Phong por facetas (CE2 de
# A. Ortiz Tello). Antes comprueba con Visibility que el objetivo está iluminado,
# en línea de visión y fuera del cono de deslumbramiento del Sol.
#
# Funciones: compute_apparent_magnitude (API), _compute_apparent_magnitude_epoch,
# _phong_model.
#
# Unidades: posiciones y distancias en km, áreas de las facetas en km² (A/d² sin
# unidades); phi_sun de photometry_config.
# =============================================================================

module Photometry

export compute_apparent_magnitude

using ..Types
using ..Constants
using ..Transformations
using ..Visibility      # eclipse, línea de visión y deslumbramiento
using ..Geometry: DebrisModel

using LinearAlgebra
using Statistics
using StaticArrays

# =============================================================================
# Public API
# =============================================================================

"""
    compute_apparent_magnitude(obs_traj, obj_traj, sun_dir, debris, phot_cfg, DCM_t) -> NamedTuple

Magnitud aparente de un objetivo visto desde un observador en cada instante.

Por qué existe: es la magnitud que la cámara convierte en flujo en la imagen y que
se guarda en las medidas LOS.

# Argumentos
- `obs_traj::TrajectoryData`, `obj_traj::TrajectoryData`: trayectorias (ECI, km).
- `sun_dir::Vector{SVector{3,Float64}}`: dirección del Sol en ECI en cada instante.
- `debris::DebrisModel`: facetas y propiedades ópticas del objetivo (`generate_debris`).
- `phot_cfg`: `photometry_config` (usa `phi_sun` y `glare_angle_deg`).
- `DCM_t::Vector{SMatrix{3,3,Float64,9}}`: DCM cuerpo -> ECI del objetivo en cada instante.

# Devuelve
- `(t, Vmag, in_eclipse, occluded, blinded, visible)`: vectores por instante;
  `Vmag = Inf` si el objetivo no se ve.

# Excepciones
- Ninguna propia; las trayectorias, `sun_dir` y `DCM_t` deben tener la misma longitud.
"""
function compute_apparent_magnitude(obs_traj::TrajectoryData, obj_traj::TrajectoryData, sun_dir::Vector{SVector{3,Float64}},
    debris::DebrisModel, phot_cfg, DCM_t::Vector{SMatrix{3,3,Float64,9}})

    N = length(obs_traj.t)

    Vmag       = Vector{Float64}(undef, N)
    in_eclipse = falses(N)
    occluded   = falses(N)
    blinded    = falses(N)
    visible    = falses(N)

    u_n_eci = similar(debris.normals)                 # normales en ECI, reutilizadas en cada instante

    @inbounds for k in 1:N

        # Rotate normals from body frame to ECI
        mul!(u_n_eci, DCM_t[k], debris.normals)

        res = _compute_apparent_magnitude_epoch(obs_traj.r[k], obj_traj.r[k], sun_dir[k],
                                                u_n_eci, debris, phot_cfg)

        Vmag[k]       = res.Vmag
        in_eclipse[k] = res.in_eclipse
        occluded[k]   = res.occluded
        blinded[k]    = res.blinded
        visible[k]    = res.visible
    end

    return (
        t = obs_traj.t,
        Vmag = Vmag,
        in_eclipse = in_eclipse,
        occluded = occluded,
        blinded = blinded,
        visible = visible
    )
end

# =============================================================================
# Internal photometry logic
# =============================================================================

"""
    _compute_apparent_magnitude_epoch(r_obs, r_obj, u_sun, u_n_eci, debris, phot_cfg) -> NamedTuple

Visibilidad y magnitud en un instante: si el objetivo está en sombra, tapado por la
Tierra o deslumbrado, `Vmag = Inf`; si no, la del modelo de Phong.

Por qué existe: separa las condiciones de visibilidad (Visibility) del modelo de
reflexión (`_phong_model`).

# Argumentos
- `r_obs`, `r_obj`: posiciones en ECI (km); `u_sun`: dirección del Sol en ECI.
- `u_n_eci::Matrix{Float64}`: normales de las facetas en ECI (3 × M).
- `debris::DebrisModel`, `phot_cfg`: como en `compute_apparent_magnitude`.

# Devuelve
- `(Vmag, in_eclipse, occluded, blinded, visible)`.

# Excepciones
- Ninguna.
"""
@inline function _compute_apparent_magnitude_epoch(r_obs::SVector{3,Float64}, r_obj::SVector{3,Float64},
    u_sun::SVector{3,Float64}, u_n_eci, debris::DebrisModel, phot_cfg)

    in_eclipse, _ = eclipseEarthECI(r_obj, u_sun)
    hasLOS = losECI(r_obs, r_obj)
    occluded = !hasLOS
    blinded, _ = sunGlareConeECI(u_sun, r_obs, r_obj, phot_cfg.glare_angle_deg)

    if in_eclipse || occluded || blinded
        return (Vmag = Inf, 
                in_eclipse = in_eclipse, 
                occluded = occluded,
                blinded = blinded, 
                visible = false)
    end

    v_obs = r_obs - r_obj
    d = norm(v_obs)
    u_obs = v_obs / d

    Vmag, has_valid_facets = _phong_model(d, u_n_eci, u_sun, u_obs, debris, phot_cfg)

    return (Vmag = Vmag, 
            in_eclipse = false, 
            occluded = false,
            blinded = false, 
            visible = has_valid_facets)
end

"""
    _phong_model(d, u_n, u_sun, u_obs, debris, phot_cfg) -> (Float64, Bool)

Flujo reflejado hacia el observador con el modelo de Phong, sumado faceta a faceta.
Para cada faceta j vista e iluminada (cos θ_in > 0 y cos θ_v > 0):

    F_j = Φ (C_diff + C_spec cos^α(θ_s) / cos θ_in) cos θ_in · A_j cos θ_v / d²

con θ_in el ángulo entre la normal n y el Sol s, θ_v entre n y el observador, y θ_s
entre el observador y la dirección especular u_spec = 2 (n·s) n - s, que ya es
unitaria: |u_spec|² = 4 (n·s)² - 4 (n·s)² + 1 = 1. La magnitud sale de F = Σ F_j.

Por qué existe: es el modelo de reflexión del objetivo; recorre las facetas en un
bucle sin arrays temporales, en el mismo orden que la ecuación.

# Argumentos
- `d::Float64`: distancia observador-objetivo (km).
- `u_n::Matrix{Float64}`: normales de las facetas en ECI (3 × M).
- `u_sun`, `u_obs`: direcciones unitarias del objetivo al Sol y al observador (ECI).
- `debris::DebrisModel`: áreas (km²), `C_spec`, `C_diff` y `alpha`.
- `phot_cfg`: `photometry_config` (usa `phi_sun`).

# Devuelve
- `(Vmag, true)`, o `(Inf, false)` si ninguna faceta está a la vez iluminada y vista.

# Excepciones
- Ninguna.
"""
function _phong_model(d::Float64, u_n, u_sun::SVector{3,Float64}, u_obs::SVector{3,Float64},
    debris::DebrisModel, phot_cfg)

    TOL = 1e-9
    F_tot = 0.0                                       # flujo total (sin dividir aún por d²)
    has_valid_facets = false

    for j in axes(u_n, 2)
        n = SVector(u_n[1, j], u_n[2, j], u_n[3, j])  # normal de la faceta j (ECI)
        cos_in = dot(n, u_sun)
        cos_v  = dot(n, u_obs)
        (cos_in > TOL && cos_v > TOL) || continue     # faceta no iluminada o no vista
        has_valid_facets = true

        u_spec = 2.0 * cos_in * n - u_sun             # dirección de reflexión especular (unitaria)
        cos_s = max(dot(u_spec, u_obs), 0.0)
        rho_spec = debris.C_spec * cos_s^debris.alpha / cos_in
        F_tot += phot_cfg.phi_sun * (debris.C_diff + rho_spec) * cos_in * debris.areas[j] * cos_v
    end
    has_valid_facets || return (Inf, false)

    F_tot /= d^2
    # [PENDIENTE] Origen de la constante -12.6988 por justificar (Vref_phong y Fref_phong no se usan)
    Vmag = -12.6988 - 2.5 * log10(F_tot)

    return (Vmag, true)
end
end # module Photometry
