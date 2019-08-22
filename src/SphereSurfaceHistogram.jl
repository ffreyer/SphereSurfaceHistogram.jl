module SphereSurfaceHistogram

using Requires

include("binning.jl")
include("utils.jl")
export SSHBinner
export partition_sphere1, partition_sphere2
export push!, append!, unsafe_append!

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
    @require GeometryTypes="4d00f742-c7ba-57c2-abde-4428a4b178cb" begin
        using .GeometryTypes
        include("extras/geom.jl")
        export face_mesh, vertex_mesh
        export bin_positions
        export line_segments_minimal, line_segments
    end
    @require Colors="5ae59095-9a9b-59fe-a467-6f913c188581" begin
        using .Colors
        include("extras/color_utils.jl")
        export to_hue, to_hue2, to_alpha, to_alpha_hue
    end
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        using .Makie
        include("extras/Makie.jl")
        export histogram, plot_debug
    end
end

end # module
