# SphereSurfaceHistogram

This package does histogram binning on a unit sphere. It can also generate
a mesh where each vertex relates to a bin.

Binning a single value takes on the order of `40ns`. 

```julia
using SphereSurfaceHistogram

# Create a Binner with approximately 1000 bins. The number of bins is optimized
# so that each bin has (approximately) the same size
N_bins = 1000
binner = SSHBinner(N_bins)

# random unit vectors are included seperately, so they aren't loaded if they
# arent needed. StaticArrays need to be installed for this.
include(SphereSurfaceHistogram.PATHS[:point_gen])
values = random_unit_vector(1000)

# values (3 component vectors) can be pushed with
push!(binner, values[1])
append!(binner, values)

# The mesh generation is sperated as well. With Makie, GeometryTypes and Colors
# is installed you can use
include(SphereSurfaceHistogram.PATHS[:Makie])

# generate a mesh where vertices correspond to bins and plot it using Makie
# the color (hue) corresponds to the relative bin filling
plot(binner)

# show the above with rectangles approximating bin areas and points marking
# the bin centers.
# note that the bins follow the curvature of the sphere, so this visualization
# is not exact.
plot_debug(binner)

# Just the bin -> vertex mesh
bin_vertex_mesh = to_dual_mesh(binner)

# just the bin -> area mesh
bin_area_mesh = to_rects(binner)

# Color mapping
colors = to_hue(binner)
```

Example plot with `1938` bins and `10 000` random unit vectors pushed. The left side shows `plot_debug(binner)`, i.e. the histogram with bin areas and centers indicated. The right side shows `plot(binner)`, i.e. colors relative to the respective bin filling.

![](https://github.com/ffreyer/SphereSurfaceHistogram.jl/blob/gh-pages/example/random_unit_vectors.png)
