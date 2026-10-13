module Transformations

export Rz, θ_g, rotate_z, llh_to_ecef, kepler_to_cartesian, theta0_from_utc, axis_angle_rotation

using ..Constants
using ..Types

using StaticArrays

@inline function θ_g(θ₀::Float64, t::Float64)
    return θ₀ + ω_earth * t
end

@inline function Rz(θ::Float64)
    c, s = cos(θ), sin(θ)
    @SMatrix [  c   s   0.0;
               -s   c   0.0;
                0.0 0.0 1.0 ]
end

@inline function rotate_z(r::SVector{3,Float64}, θ::Float64)
    return transpose(Rz(θ)) * r
end

@inline function llh_to_ecef(lat_rad::Float64, lon_rad::Float64, h_km::Float64)
    sφ, cφ = sin(lat_rad), cos(lat_rad)
    N = WGS84_A / sqrt(1 - WGS84_E2*sφ^2)
    x = (N + h_km) * cφ * cos(lon_rad)
    y = (N + h_km) * cφ * sin(lon_rad)
    z = (N*(1 - WGS84_E2) + h_km) * sφ
    return SVector{3,Float64}(x, y, z)
end


@inline function axis_angle_rotation(axis::SVector{3,Float64}, angle::Float64)::SMatrix{3,3,Float64,9}

    ux, uy, uz = axis

    c = cos(angle)
    s = sin(angle)
    d = 1.0 - c

    return @SMatrix [
        c + ux*ux*d       ux*uy*d - uz*s   ux*uz*d + uy*s;
        uy*ux*d + uz*s    c + uy*uy*d      uy*uz*d - ux*s;
        uz*ux*d - uy*s    uz*uy*d + ux*s   c + uz*uz*d
    ]
end

#------------------------------------FUNCTIONS----------------------------------------------

function kepler_to_cartesian(dict)::State

    # Extraer campos EXACTOS del JSON
    a = Float64(dict.SMA)
    e = Float64(dict.ECC)
    i = deg2rad(Float64(dict.INC))
    Ω = deg2rad(Float64(dict.RAAN))
    ω = deg2rad(Float64(dict.AOP))
    ν = deg2rad(Float64(dict.TA))

    p = a * (1 - e^2)

    # Perifocal
    r_pf = @SVector [
        p * cos(ν) / (1 + e*cos(ν)),
        p * sin(ν) / (1 + e*cos(ν)),
        0.0
    ]

    v_pf = @SVector [
        -sqrt(μ/p) * sin(ν),
         sqrt(μ/p) * (e + cos(ν)),
         0.0
    ]

    # Rotaciones
    R3Ω = @SMatrix [
        cos(Ω) -sin(Ω) 0.0;
        sin(Ω)  cos(Ω) 0.0;
        0.0     0.0    1.0
    ]

    R1i = @SMatrix [
        1.0 0.0      0.0;
        0.0 cos(i) -sin(i);
        0.0 sin(i)  cos(i)
    ]

    R3ω = @SMatrix [
        cos(ω) -sin(ω) 0.0;
        sin(ω)  cos(ω) 0.0;
        0.0     0.0    1.0
    ]

    Q = R3Ω * R1i * R3ω

    return State(Q * r_pf, Q * v_pf)
end


function theta0_from_utc(jd::Float64)
    era = 2π*(0.7790572732640 + 1.00273781191135448*(jd - jd_J2000))
    return mod(era, 2π)
end
end