
struct SSHBinner <: AbstractSSH
    tessellation::SphereTessellationMap
    bins::Vector{Int64}
end

function SSHBinner(tess::SphereTessellationMap)
    return SSHBinner(tess, zeros(Int64, length(tess)))
end

"""
    SSHBinner(N_bins::Integer)

Generate a Sphere Surface Histogram Binner with approximately N_bins.
"""
SSHBinner(Nbins; kwargs...) = SSHBinner(SphereTessellationMap(Nbins; kwargs...))

"""
    SSHBinner(bins)

Generate a Sphere Surface Histogram Binner using the `bins` of another
SSHBinner. This effectively creates a copy of the other source `SSHBinner`.
"""
function SSHBinner(bins::Vector; method=partition_sphere2)
    N_bins = length(bins)
    thetas, phi_divisions = method(4pi/N_bins)
    tess = SphereTessellationMap(thetas, phi_divisions)
    return SSHBinner(tess, bins)
end

"""
    push!(B::SSHBinner, vec)

Bins a single cartesian vec (three dimensional unit vector).
"""
function Base.push!(B::SSHBinner, vec::AbstractVector)
    B.bins[bin_index(B.tessellation, vec)] += 1
    return nothing
end

"""
    push!(B::SSHBinner, theta, phi)

Bins a single unit vector described by two angles theta and phi. The angles
follow the physics convention, i.e.
- `theta` is 0 in +z direction and pi at -z direction and
- `phi` is 0 at +x direction, increasing counterclockwise (pi/2 at +y)
"""
function Base.push!(B::SSHBinner, theta::Real, phi::Real)
    B.bins[bin_index(B.tessellation, theta, phi)] += 1
    return nothing
end

get_value(B::SSHBinner, idx::Int) = B.bins[idx]
get_values(B::SSHBinner) = B.bins
