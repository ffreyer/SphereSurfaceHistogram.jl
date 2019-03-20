module SphereSurfaceHistogram

include("utils.jl")
# export binary_find, to_cartesian

include("binning.jl")
export SSHBinner
export bin!, bin_many!

#include("Makie.jl")
_DIR = @__DIR__
const PATHS = Dict(
    :point_gen => _DIR * "/extras/point_gen.jl",
    :Makie => _DIR * "/extras/Makie.jl"
)

end # module
