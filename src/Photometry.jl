module Photometry

export compute_apparent_magnitude

using ..Types
using ..Constants
using ..Transformations
using ..Environment

using LinearAlgebra
using Statistics
using StaticArrays

const DCM_I = @SMatrix [
    1.0 0.0 0.0;
    0.0 1.0 0.0;
    0.0 0.0 1.0
]

# =============================================================================
# Public API
# =============================================================================

function compute_apparent_magnitude(obs_traj::TrajectoryData, obj_traj::TrajectoryData, sun_dir::Vector{SVector{3,Float64}},
    debris::Dict{Symbol,Any}, phot_cfg, DCM_t::Vector{SMatrix{3,3,Float64,9}})

    N = length(obs_traj.t)

    Vmag       = Vector{Float64}(undef, N)
    in_eclipse = falses(N)
    occluded   = falses(N)
    blinded    = falses(N)
    visible    = falses(N)

    u_n_body = debris[:normals]
    u_n_eci  = similar(u_n_body)

    A      = debris[:areas]
    alpha  = debris[:alpha]
    C_spec = debris[:C_spec]
    C_diff = debris[:C_diff]

    @inbounds for k in 1:N

        # Rotate normals from body frame to ECI
        mul!(u_n_eci, DCM_t[k], u_n_body)

        res = _compute_apparent_magnitude_epoch(
            obs_traj.r[k],
            obj_traj.r[k],
            sun_dir[k],
            u_n_eci,
            A,
            alpha,
            C_spec,
            C_diff,
            phot_cfg
        )

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

@inline function _compute_apparent_magnitude_epoch(r_obs::SVector{3,Float64}, r_obj::SVector{3,Float64},
    u_sun::SVector{3,Float64}, u_n_eci, A, alpha::Float64, C_spec::Float64, C_diff::Float64,
    phot_cfg)

    in_eclipse, _, _ = eclipseEarthECI(r_obj, u_sun)
    hasLOS, _ = losECI(r_obs, r_obj)
    occluded = !hasLOS
    blinded, _, _ = sunGlareConeECI(u_sun, r_obs, r_obj, phot_cfg.glare_angle_deg)

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

    Vmag, has_valid_facets = _phong_model(d, C_spec, C_diff, alpha, u_n_eci, u_sun, u_obs, A, phot_cfg)

    return (Vmag = Vmag, 
            in_eclipse = false, 
            occluded = false,
            blinded = false, 
            visible = has_valid_facets)
end

function _phong_model(d::Float64, C_spec::Float64, C_diff::Float64, alpha::Float64,
    u_n, u_sun::SVector{3,Float64}, u_obs::SVector{3,Float64}, A, phot_cfg)

    TOL = 1e-9

    cos_in = vec(sum(u_n .* u_sun, dims = 1))
    cos_v  = vec(sum(u_n .* u_obs, dims = 1))

    ok = (cos_in .> TOL) .& (cos_v .> TOL)
    any(ok) || return (Inf, false)

    u_n_ok    = u_n[:, ok]
    cos_in_ok = cos_in[ok]
    cos_v_ok  = cos_v[ok]
    A_ok      = A[ok]

    u_spec = 2.0 .* (u_n_ok .* reshape(cos_in_ok, 1, :)) .- reshape(u_sun, 3, 1)
    u_spec ./= sqrt.(sum(abs2, u_spec; dims = 1))

    cos_s = max.(vec(sum(u_spec .* u_obs, dims = 1)), 0.0)

    rho_spec = C_spec .* (cos_s .^ alpha) ./ cos_in_ok
    F_sun    = phot_cfg.phi_sun .* (C_diff .+ rho_spec) .* cos_in_ok
    F_obs    = F_sun .* A_ok .* cos_v_ok ./ (d^2)

    F_tot = sum(F_obs)
    # Vmag = -26.7 - 2.5 * log10(F_tot / phi)
    Vmag = -12.6988 - 2.5 * log10(F_tot)

    return (Vmag, true)
end
end # module Photometry