# SphereSurfaceHistogram.jl

`SphereSurfaceHistogram` is a package for binning 3D unit vectors.
Each bin has the same size - i.e. it covers (approximately) the same (curved) area on a unit sphere.

There are currently two ways to bin values in the package.
The first is `SSHBinner` which simply counts how often a bin has been pushed to.
The other is `SSHAverager` which reports the average float value that has been pushed to a specific bin.

```@contents
```

## Quickstart


An empty binner with approximately `10_000` bins can be created with

```@repl
using SphereSurfaceHistogram
binner = SSHBinner(10_000)
```

Adding vectors can be done using `push!` or `append!`

```@repl
using SphereSurfaceHistogram
binner = SSHBinner(10_000)
push!(binner, [1.0, 0.0, 0.0])
append!(binner, [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
```

!!! warning
    SphereSurfaceHistogram does not verify that the given vector is normalized.

`SSHAverager` works analogously:

```@repl
using SphereSurfaceHistogram
binner = SSHaverager(10_000)
push!(binner, [1.0, 0.0, 0.0], 1.6)
append!(binner, [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]], [0.3, 0.6])
```

You can also push to a binner using angles of a spherical coordinate system.
This follows the physics convention as mentioned on [wikipedia](https://en.wikipedia.org/wiki/Spherical_coordinate_system), i.e. with $theta \in (0, pi)$ being the angle between the z axis and the vector and $phi \in (0, 2pi)$ extending counterclockwise from the x axis.

```@repl
push!(binner, 0.5pi, 1.7pi) # theta, phi
```

!!! warning
    SphereSurfaceHistogram does not verify that the given angles are in bounds.

## Accessing the Histogram

SphereSurfaceHistogram provides custom getindex methods working with angles in radians.

```@docs
getindex(::SSHBinner, ::Real, ::Real)
```

## Plotting

SphereSurfaceHistogram includes a couple of Makie recipes. Most notably is the `histogram` recipe. It creates a mesh where each vertex represents a bin. Colors are then picked based on the counts in each bin.

![Example Plot](assets/main_example.png)

The image above uses the following code

```julia
using SphereSurfaceHistogram, GLMakie
binner1 = SSHBinner(10_000)
binner2 = SSHBinner(10_000)
for _ in 1:1000
    append!(binner1, random_unit_vector(1_000))
end
for _ in 1:1_000_000
    push!(binner2, normalize(2rand(3) .- 1.0))
end
fig = Figure()
histogram(fig[1, 1], binner1)
histogram(fig[1, 2], binner2)
fig
```

!!! info
    Plotting methods are not loaded by default. SphereSurfaceHistogram makes use of `Requires.jl` to dynamically load plotting related methods when `Makie.jl` becomes available.

!!! info
    Mesh generation methods can also be loaded individually. If you want to use the mesh generation methods without `Makie`, you should import `GeometryBasics.jl`.