using GeometryTypes, Colors, Makie
include("Makie_utils.jl")
include("Makie_geom.jl")

"""
    histogram(B::SSHBinner[, color_method = to_hue])

Plots a sphere where each bin filling is mapped to a color through color_method.
"""
function histogram(B::SSHBinner, color_method=to_hue)
    _colors = color_method(B)
    Makie.mesh(to_dual_mesh(B), color = _colors, transparency = false)
end


"""
    plot_debug(B::SSHBinner[, color_method = to_hue])

Plots a sphere with colors representing bin fillings, approximate bin area
(rectangles) and points at the bin centers.
"""
function plot_debug(B::SSHBinner, color_method = to_hue)
    _colors = color_method(B)
    scene = Makie.mesh(to_dual_mesh(B), color = _colors, transparency = false)
    Makie.linesegments!(scene, 1.001line_segments(B), color = :black)
    # Makie.mesh!(scene1, to_rects(B, extrude), color = RGBA(1, 1, 1, 0.2), transparency = true)
    Makie.scatter!(scene, 1.005dual_points(B), color=:black, markersize=0.01)
    scene
end
