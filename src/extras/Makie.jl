"""
    histogram(B::SSHBinner[, color_method = to_hue])

Plots a sphere where each bin filling is mapped to a color through color_method.
"""
function histogram(B::SSHBinner, color_method=to_hue)
    _colors = color_method(B)
    is_transparent = any(c -> c.alpha < 1.0, _colors)
    Makie.mesh(vertex_mesh(B), color = _colors, transparency = is_transparent)
end


"""
    plot_debug(B::SSHBinner[, color_method = to_hue])

Plots a sphere with colors representing bin fillings, approximate bin area
(rectangles) and points at the bin centers.
"""
function plot_debug(B::SSHBinner, color_method = to_hue)
    _colors = color_method(B)
    is_transparent = any(c -> c.alpha < 1.0, _colors)
    scene = Makie.mesh(vertex_mesh(B), color = _colors, transparency = is_transparent)
    Makie.linesegments!(scene, 1.001line_segments(B), color = :black)
    # Makie.mesh!(scene1, to_rects(B, extrude), color = RGBA(1, 1, 1, 0.2), transparency = true)
    Makie.scatter!(scene, 1.005bin_positions(B), color=:black, markersize=0.01)
    scene
end
