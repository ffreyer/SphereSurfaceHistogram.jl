module SphereSurfaceHistogram

using Reexport
@reexport using StaticArrays
import Statistics: mean, std

include("utils.jl")
export rand_spin, binary_find, to_cartesian

include("binning.jl")
export SSHBinner
export bin!, bin_many!

include("Makie.jl")

end # module
