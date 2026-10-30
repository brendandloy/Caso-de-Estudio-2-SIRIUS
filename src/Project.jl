module Project

# ========================
# submodul inclussions
# ========================

include("Types.jl")
include("Constants.jl")
include("Parameters.jl")          # parámetros maestros del programa (no los elige el usuario)
include("Transformations.jl")
include("Attitude/Attitude.jl")     # actitud de los objetivos; Environment la usa al leer el JSON
include("Scenario/Scenario.jl")     # generador spec -> JSON completo (lo usa scripts/main.jl)
include("Environment.jl")
include("Propagator.jl")
include("Visualization.jl")
include("Camera_sim.jl")
include("Geometry.jl")
include("Photometry.jl")

# ========================
# exported submodules
# ========================

using .Types
using .Constants
using .Parameters
using .Transformations
using .Attitude
using .Scenario
using .Environment
using .Propagator
using .Visualization
using .CameraSim
using .Geometry
using .Photometry

end # module Project