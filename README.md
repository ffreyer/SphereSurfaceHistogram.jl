# SphereSurfaceHistogram

[![](https://github.com/ffreyer/SphereSurfaceHistogram.jl/workflows/Documentation/badge.svg)](https://ffreyer.github.io/SphereSurfaceHistogram.jl/dev)

This package does histogram binning on a unit sphere. It can also generate
a mesh where each vertex relates to a bin.

Binning a single value takes on the order of `30ns`.

---

Example plot with `1938` bins and `10 000` random unit vectors pushed. The left side shows `plot_debug(binner)`, i.e. the histogram with bin areas and centers indicated. The right side shows `histogram(binner)`, i.e. colors relative to the respective bin filling.

![](https://github.com/ffreyer/SphereSurfaceHistogram.jl/blob/gh-pages/example/random_unit_vectors.png)

---

```julia
using SphereSurfaceHistogram

# random_unit_vector() is an utility function that is only loaded when
# StaticArrays is loaded. We'll use it to gnerate some sample data
using StaticArrays
values = random_unit_vector(1000)



# Create a Binner with approximately 1000 bins. The number of bins is optimized
# so that each bin has (approximately) the same size
N_bins = 1000
binner = SSHBinner(N_bins)

# values can be pushed (binned) with
push!(binner, values[1])    # one 3-component vector
append!(binner, values)     # many 3-component vectors

# Bins can be accessed via indexing B[theta, phi]
B[pi/2, pi]           # bin closest to (-1, 0, 0)
B[pi/2, :]            # phi-ring closest to z = 0
B[:, pi/2]            # theta-ring closest to x = 0, y = 1
B[:, :]               # all bins, indexed as [theta_index][phi_index]
B[pi/2, 0:pi/2:2pi]   # bins closest to (cos(phi), sin(phi), 0) for phi = 0:pi/2:2pi


# Additionally SSHBinner can be supplied with a `method`. There are two options,
# `partition_sphere1` and `partition_sphere2`. Both partition the sphere in theta
# first, creating some number of rings with phi = 0..2pi. The former then creates
# any number bins for each ring, the latter is restricted to some power of two,
# which means that bins a re more frequently aligned.
# `partition_sphere2` is the default.
B = SSHBinner(N_bins, method=partition_sphere1)
```

---

Additionally, this package includes some plotting utilities. Just like `random_unit_vector`, they are only loaded when the appropriate packages are
loaded. Continueing the example above:

```julia
using Makie


# Plot a histogram
# For this a mesh is generated where each vertex corresponds to one bin. The
# vertices are then colored based on the correspondign bin value.
histogram(binner, colormap=:matter, :colorscale=:relative, outline=:false)
```

---

Some other exportet functions (from plotting utilities) include:

```julia
# Generate a mesh where each vertex corresponds to the center of a bin
bin_vertex_mesh = vertex_mesh(binner)

# Generate a disconnected mesh where each face corresponds to a bin
bin_area_mesh = face_mesh(binner)

# smoothed line segments, outlining each bin
_lines = line_segments(binner)
# unsmoothed version
_lines = line_segments_minimal(binner)

# positions of each bin center
points = bin_positions(binner)
```
