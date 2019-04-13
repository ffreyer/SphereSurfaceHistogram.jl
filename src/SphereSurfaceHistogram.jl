module SphereSurfaceHistogram

include("binning.jl")
export SSHBinner
export partition_sphere1, partition_sphere2
export push!, append!, unsafe_append!

#include("Makie.jl")
_DIR = @__DIR__
const PATHS = Dict(
    :point_gen => _DIR * "/extras/point_generator.jl",
    :Makie => _DIR * "/extras/Makie.jl"
)

end # module
