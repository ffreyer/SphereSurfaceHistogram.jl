module SphereSurfaceHistogram

include("utils.jl")
# export binary_find, to_cartesian

include("binning.jl")
export SSHBinner
export push!, append!

#include("Makie.jl")
_DIR = @__DIR__
const PATHS = Dict(
    :point_gen => _DIR * "/extras/point_generator.jl",
    :Makie => _DIR * "/extras/Makie.jl"
)

end # module
