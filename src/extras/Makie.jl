import .Makie: convert_arguments, plot!

# bin centers
function convert_arguments(P::Type{<:Scatter}, B::SSHBinner)
    convert_arguments(P, bin_positions(B))
end
# mesh with vertices @ bin centers (makes color mapping easy)
function convert_arguments(P::Type{<:Makie.Mesh}, B::SSHBinner)
    convert_arguments(P, vertex_mesh(B))
end
# boundaries of bins (slightly enlarged to work with mesh)
function convert_arguments(P::Type{<:LineSegments}, B::SSHBinner)
    convert_arguments(P, 1.001line_segments(B))
end


# Not sure if I wanna keep this. Would probably be nice to have a colorlegend
# by default here
"""
    histogram(B::SSHBinner)
    histogram!(scene, B::SSHBinner)

Plots a histogram for the given SSHBinner.

# Attributes

- `outline = false`: Draws lines around each bin
- `colormap = :matter`: Colormap used to visualize bins
- `colorscale = :absolute`: Scale used when mapping bins to colors. Can be
:relative (sets colorrange = extrema(B.bins)), :absolute (sets only the max) or
:static (colorrange remains unchanged)
- `colorrange = [0, maxiumum(b.bins)]`: Range of bin fillings represented by the
colormap.
- `transparency = false`: Set this to true if your colormap is tranparent
- `linecolor = :black`: Color of outline.
"""
@recipe(Histogram, B) do scene
    Theme(
        outline = false,
        colormap = :matter,
        transparency = false,
        linecolor = :black,
        colorscale = :absolute,
        colorrange = (0, 1)
    )
end

function map_bin_to_color(binner, colormap, colorscale, colorrange)
    map(binner, colormap, colorscale) do binner, cm, scale
        ex = extrema(binner.bins)
        if scale == :absolute
            colorrange[] = (0, ex[2])
        elseif scale == :static
        else
            colorrange[] = (ex[1], ex[2])
        end
        map(binner.bins) do b
            RGBAf0(Makie.interpolated_getindex(
                to_colormap(cm), Float64(b), colorrange[]
            ))
        end
    end
end

function plot!(plot::Histogram)
    colors = map_bin_to_color(
        plot[:B], plot[:colormap], plot[:colorscale], plot[:colorrange]
    )

    mesh!(plot, plot[:B], color = colors, transparency = plot[:transparency])
    linesegments!(plot, plot[:B], color = plot[:linecolor], visible=plot[:outline])

    plot
end


# """
#     histogram(B::SSHBinner[, color_method = to_hue])
#
# Plots a sphere where each bin filling is mapped to a color through color_method.
# """
# function histogram(B::SSHBinner, color_method=to_hue)
#     _colors = color_method(B)
#     is_transparent = any(c -> c.alpha < 1.0, _colors)
#     AbstractPlotting.mesh(vertex_mesh(B), color = _colors, transparency = is_transparent)
# end
#
#
# """
#     plot_debug(B::SSHBinner[, color_method = to_hue])
#
# Plots a sphere with colors representing bin fillings, approximate bin area
# (rectangles) and points at the bin centers.
# """
# function plot_debug(B::SSHBinner, color_method = to_hue)
    # _colors = color_method(B)
    # is_transparent = any(c -> c.alpha < 1.0, _colors)
    # scene = AbstractPlotting.mesh(vertex_mesh(B), color = _colors, transparency = is_transparent)
    # AbstractPlotting.linesegments!(scene, 1.001line_segments(B), color = :black)
    # # AbstractPlotting.mesh!(scene1, to_rects(B, extrude), color = RGBA(1, 1, 1, 0.2), transparency = true)
    # AbstractPlotting.scatter!(scene, 1.005bin_positions(B), color=:black, markersize=0.01)
    # scene
# end
