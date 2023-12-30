
struct SSHBinner <: AbstractSSH
    thetas::Vector{Float64}
    phi_divisions::Vector{Int64}

    # bin_many helpers
    zs::Vector{Float64}
    phi_N_over_2pi::Vector{Float64}
    phi_cumsum::Vector{Int64}
    N_thetas::Float64
    theta_factor::Float64

    bins::Vector{Int64}
end

"""
    SSHBinner(N_bins::Integer)

Generate a Sphere Surface Histogram Binner with approximately N_bins.
"""
function SSHBinner() end

"""
    SSHBinner(bins)

Generate a Sphere Surface Histogram Binner using the `bins` of another
SSHBinner. This effectively creates a copy of the other source `SSHBinner`.
"""
function SSHBinner(bins::Vector; method=partition_sphere2)
    N_bins = length(bins)
    thetas, phi_divisions = method(4pi/N_bins)

    PI_PLUS = 2.0pi + 1e-14
    phi_N_over_2pi = phi_divisions ./ PI_PLUS
    zs = cos.(thetas)
    zs[end] -= 1e-14
    phi_cumsum = cumsum(phi_divisions)
    N_thetas = Float64(length(thetas)-1)
    theta_factor = Float64(length(thetas)-1) / (pi + 1e-14)

    return SSHBinner(thetas, phi_divisions, zs, phi_N_over_2pi, phi_cumsum, N_thetas, theta_factor, bins)
end

"""
    push!(B::SSHBinner, vec)

Bins a single cartesian vec (three dimensional unit vector).
"""
function Base.push!(B::SSHBinner, vec::AbstractVector)
    B.bins[bin_index(B, vec)] += 1
    return nothing
end

"""
    push!(B::SSHBinner, theta, phi)

Bins a single unit vector described by two angles theta and phi. The angles
follow the physics convention, i.e.
- `theta` is 0 in +z direction and pi at -z direction and
- `phi` is 0 at +x direction, increasing counterclockwise (pi/2 at +y)
"""
function Base.push!(B::AbstractSSH, theta::AbstractFloat, phi::AbstractFloat)
    B.count[bin_index(B, theta, phi)] += 1
    return nothing
end

get_value(B::SSHBinner, idx::Int) = B.bins[idx]
get_values(B::SSHBinner) = B.bins
