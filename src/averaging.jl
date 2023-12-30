
struct SSHAverager <: AbstractSSH
    thetas::Vector{Float64}
    phi_divisions::Vector{Int64}

    # bin_many helpers
    zs::Vector{Float64}
    phi_N_over_2pi::Vector{Float64}
    phi_cumsum::Vector{Int64}
    N_thetas::Float64
    theta_factor::Float64

    sums::Vector{Float64}
    counts::Vector{Int64}
end

function SSHAverager(thetas, phi_divisions, zs, phi_N_over_2pi, phi_cumsum, N_thetas, theta_factor, counts)
    sums = zeros(Float64, length(counts))
    return SSHAverager(thetas, phi_divisions, zs, phi_N_over_2pi, phi_cumsum, N_thetas, theta_factor, sums, counts)
end

"""
    SSHAverager(N_bins::Integer)

Generate a Sphere Surface Histogram Averager with approximately N_bins.
"""
function SSHAverager() end

"""
    SSHAverager(sums, counts)

Generate a Sphere Surface Histogram Averager using the `bins` of another
SSHAverager. This effectively creates a copy of the other source `SSHAverager`.
"""
function SSHAverager(sums::Vector{Float64}, counts::Vector{Int64}; method=partition_sphere2)
    N_bins = length(sums)
    thetas, phi_divisions = method(4pi/N_bins)

    PI_PLUS = 2.0pi + 1e-14
    phi_N_over_2pi = phi_divisions ./ PI_PLUS
    zs = cos.(thetas)
    zs[end] -= 1e-14
    phi_cumsum = cumsum(phi_divisions)
    N_thetas = Float64(length(thetas)-1)
    theta_factor = Float64(length(thetas)-1) / (pi + 1e-14)

    return SSHAverager(thetas, phi_divisions, zs, phi_N_over_2pi, phi_cumsum, N_thetas, theta_factor, sums, counts)
end

"""
    push!(B::SSHAverager, vec::AbstractVector, value::AbstractFloat)

Adds a single value at a position given by a cartesian vec (three dimensional
unit vector) to the mean.
"""
function Base.push!(B::SSHAverager, vec::AbstractVector, value::AbstractFloat)
    idx = bin_index(B, vec)
    B.sums[idx] += value
    B.counts[] += 1
    return nothing
end

"""
    push!(B::SSHAverager, theta::AbstractFloat, phi::AbstractFloat, value::AbstractFloat)

Adds a single value at the given angles theta and phi to the running mean. The
angles follow the physics convention, i.e.
- `theta` is 0 in +z direction and pi at -z direction and
- `phi` is 0 at +x direction, increasing counterclockwise (pi/2 at +y)
"""
function Base.push!(B::AbstractSSH, theta::AbstractFloat, phi::AbstractFloat, value::AbstractFloat)
    idx = bin_index(B, theta, phi)
    B.sums[idx] += value
    B.counts[idx] += 1
    return nothing
end

get_value(B::SSHAverager, idx::Int) = B.sums[idx] / B.counts[idx]
get_values(B::SSHAverager) = B.sums ./ B.counts