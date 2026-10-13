module Visualization

export setup_scene, add_orbit!

using ..Constants

using GLMakie
using StaticArrays
using GeometryBasics: Point3

#------------------------------------FUNCTIONS----------------------------------------------

function setup_scene()
    set_theme!(theme_dark())
    fig = Figure(size=(800, 800))
    ax = Axis3(fig[1, 1], 
        xlabel="x [10³ km]", 
        ylabel="y [10³ km]", 
        zlabel="z [10³ km]")
    ax.aspect = :data

    # Tierra (esfera)
    R_scaled = R / 1000
    θ = range(0, 2π, length=50)
    φ = range(0, π, length=50)
    x = [R_scaled * sin(ϕ) * cos(θi) for ϕ in φ, θi in θ]
    y = [R_scaled * sin(ϕ) * sin(θi) for ϕ in φ, θi in θ]
    z = [R_scaled * cos(ϕ) for ϕ in φ, θi in θ]
    surface!(ax, x, y, z, colormap=:blues, transparency=true, alpha=0.7)

    # Ecuador (sin transparencia)
    l = R_scaled * 1.5
    ecuador = [
        Point3f(-l, -l, 0),
        Point3f(l, -l, 0),
        Point3f(l, l, 0),
        Point3f(-l, l, 0)
    ]
    poly!(ax, ecuador, color=:red, transparency=true, alpha=0.4)

    return fig, ax
end


function add_orbit!(ax, traj::Vector{SVector{3, Float64}}, label::String)
    pts = [Point3f(s./ 1000) for s in traj]  # Convertir km a 10³ km
    lines!(ax, pts, label=label, linewidth=2)
end

end