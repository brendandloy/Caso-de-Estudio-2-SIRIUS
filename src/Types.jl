module Types

export State, TrajectoryData, GroundStationLLH, StationTrack, LOSMeasurement, AttitudeData

using StaticArrays

struct State
    r::SVector{3, Float64}
    v::SVector{3, Float64}
end

struct TrajectoryData
    r::Vector{SVector{3, Float64}}  # posiciones
    v::Vector{SVector{3, Float64}}  # velocidades
    t::Vector{Float64}              # tiempos
end

struct GroundStationLLH
    lat_deg::Float64
    lon_deg::Float64
    h_m::Float64               # altura elipsoidal (m)
end

struct StationTrack
    r::Vector{SVector{3,Float64}}  # posiciones ECI
    t::Vector{Float64}             # tiempos (s)
end

struct LOSMeasurement
    t::Float64
    dir_true::SVector{3,Float64}
    dir_meas::SVector{3,Float64}
    Vmag::Float64
end

struct AttitudeData
    DCM_b2i::Vector{SMatrix{3,3,Float64,9}}
end

end