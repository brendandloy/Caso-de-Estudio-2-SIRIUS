module Project

# ========================
# submodul inclussions
# ========================

include("Types.jl")
include("Constants.jl")
include("Transformations.jl")
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
using .Transformations
using .Environment
using .Propagator
using .Visualization
using .CameraSim
using .Geometry
using .Photometry

end # module Project