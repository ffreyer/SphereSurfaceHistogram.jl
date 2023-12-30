
struct SSHAverager <: AbstractSSH
    tessellation::SphereTessellationMap
    sums::Vector{Float64}
    counts::Vector{Int64}
end

function SSHAverager(tess::SphereTessellationMap)
    return SSHAverager(tess, zeros(Float64, length(tess)), zeros(Int64, length(tess)))
end

"""
    SSHAverager(N_bins::Integer)

Generate a Sphere Surface Histogram Averager with approximately N_bins.
"""
SSHAverager(Nbins; kwargs...) = SSHAverager(SphereTessellationMap(Nbins; kwargs...))

"""
    SSHAverager(sums, counts)

Generate a Sphere Surface Histogram Averager using the `bins` of another
SSHAverager. This effectively creates a copy of the other source `SSHAverager`.
"""
function SSHAverager(sums::Vector{Float64}, counts::Vector{Int64}; method=partition_sphere2)
    N_bins = length(sums)
    thetas, phi_divisions = method(4pi/N_bins)
    tess = SphereTessellationMap(thetas, phi_divisions)
    return SSHAverager(tess, sums, counts)
end

"""
    push!(B::SSHAverager, vec::AbstractVector, value::Real)

Adds a single value at a position given by a cartesian vec (three dimensional
unit vector) to the mean.
"""
function Base.push!(B::SSHAverager, vec::AbstractVector, value::Real)
    idx = bin_index(B.tessellation, vec)
    B.sums[idx] += value
    B.counts[idx] += 1
    return nothing
end

"""
    push!(B::SSHAverager, theta::Real, phi::Real, value::Real)

Adds a single value at the given angles theta and phi to the running mean. The
angles follow the physics convention, i.e.
- `theta` is 0 in +z direction and pi at -z direction and
- `phi` is 0 at +x direction, increasing counterclockwise (pi/2 at +y)
"""
function Base.push!(B::SSHAverager, theta::Real, phi::Real, value::Real)
    idx = bin_index(B.tessellation, theta, phi)
    B.sums[idx] += value
    B.counts[idx] += 1
    return nothing
end

get_value(B::SSHAverager, idx::Int) = B.sums[idx] / max(1, B.counts[idx])
get_values(B::SSHAverager) = B.sums ./ max.(1, B.counts)