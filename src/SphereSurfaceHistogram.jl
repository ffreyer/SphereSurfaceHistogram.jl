module SphereSurfaceHistogram

using Reexport
@reexport using StaticArrays
import Statistics: mean, std

include("utils.jl")
export rand_spin, binary_find, to_cartesian

include("binning.jl")
export partition_sphere2, partition_sphere_optim2
export bin_many!

include("Makie.jl")

end # module
