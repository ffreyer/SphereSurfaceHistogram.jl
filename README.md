# SphereSurfaceHistogram

This package does histogram binning on a unit sphere. It can also generate
a mesh where each vertex relates to a bin.

Binning a single value takes on the order of `30ns`. 

---

Example plot with `1938` bins and `10 000` random unit vectors pushed. The left side shows `plot_debug(binner)`, i.e. the histogram with bin areas and centers indicated. The right side shows `histogram(binner)`, i.e. colors relative to the respective bin filling.

![](https://github.com/ffreyer/SphereSurfaceHistogram.jl/blob/gh-pages/example/random_unit_vectors.png)

---

```julia
using SphereSurfaceHistogram

# random unit vectors are included seperately, so they aren't loaded if they
# arent needed. StaticArrays need to be installed for this.
include(SphereSurfaceHistogram.PATHS[:point_gen])
values = random_unit_vector(1000)



# Binning

# Create a Binner with approximately 1000 bins. The number of bins is optimized
# so that each bin has (approximately) the same size
N_bins = 1000
binner = SSHBinner(N_bins)

# values can be pushed (binned) with
push!(binner, values[1])    # one 3-component vector
append!(binner, values)     # many 3-component vectors

# Bins can be accessed via indexing B[phi, theta]
B[pi, pi/2]           # bin closest to (-1, 0, 0)
B[:, pi/2]            # phi-ring closest to z = 0
B[pi/2, :]            # theta-ring closest to x = 0, y = 1
B[:, :]               # all bins, indexed as [theta_index][phi_index]
B[0:pi/2:2pi, pi/2]   # bins closest to (cos(phi), sin(phi), 0) for phi = 0:pi/2:2pi


# Additionally SSHBinner can be supplied with a `method`. There are two options,
# `partition_sphere1` and `partition_sphere2`. Both partition the sphere in theta
# first, creating some number of rings with phi = 0..2pi. The former then creates
# any number bins for each ring, the latter is restricted to some power of two,
# which means that bins a re more frequently aligned.
# `partition_sphere2` is the default.
B = SSHBinner(N_bins, method=partition_sphere1)



# Plotting

# The mesh generation is sperated as well. With Makie, GeometryTypes and Colors
# is installed you can use
include(SphereSurfaceHistogram.PATHS[:Makie])

# generate a mesh where vertices correspond to bins and plot it using Makie
# the color (hue) corresponds to the relative bin filling
histogram(binner)

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
