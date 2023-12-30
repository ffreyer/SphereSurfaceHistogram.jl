module SphereSurfaceHistogram

using Requires

include("generic.jl")
include("counting.jl")
include("averaging.jl")
include("utils.jl")

export SSHBinner, SSHAverager
export partition_sphere1, partition_sphere2
export unsafe_append!
export get_values

#include("Makie.jl")
# _DIR = @__DIR__
# const PATHS = Dict(
#     :point_gen => _DIR * "/extras/point_generator.jl",
#     :Makie => _DIR * "/extras/Makie.jl"
# )

function __init__()
    @require StaticArrays="90137ffa-7385-5640-81b9-e52037218182" begin
        using .StaticArrays
        include("extras/point_generator.jl")
        export random_unit_vector
    end
    @require GeometryBasics="5c1252a2-5f33-56bf-86c9-59e7332b4326" begin
        using .GeometryBasics
        include("extras/geom.jl")
        export face_mesh, vertex_mesh
        export bin_positions
        export line_segments_minimal, line_segments
    end
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        using .Makie
        include("extras/Makie.jl")
        export histogram, histogram!
    end
end

end # module
