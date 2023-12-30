################################################################################
### SSHBinner utilities (things from Base and more)
################################################################################


function Base.show(io::IO, B::T) where {T <: AbstractSSH}
    print(io, "$T with $(length(B)) bins")
end

Base.summary(io::IO, B::AbstractSSH) = print(io, "$(typeof(B)){$(length(B))}")

# "Return the number of elements in the collection."
# would be sum(B.bins)
# Maybe length(B.bins) makes more sense?
Base.length(B::SSHBinner) = length(B.bins)
Base.length(B::SSHAverager) = length(B.sums)

# TODO: What should this be exactly?
# - length per theta? -> phi_divisions
# - summarized length per theta? (1x 1, 4x 4, 20x 32, ...)
Base.size(B::AbstractSSH) = B.phi_divisions

# Clears bins
function Base.empty!(B::SSHBinner)
    B.bins .= 0
    return B
end
function Base.empty!(B::SSHAverager)
    B.sums .= 0.0
    B.counts .= 0
    return B
end

Base.maximum(B::SSHBinner) = maximum(B.bins)
Base.minimum(B::SSHBinner) = minimum(B.bins)

Base.maximum(B::SSHAverager) = maximum(B.sums)
Base.minimum(B::SSHAverager) = minimum(B.sums)


################################################################################
### getindex(B, phi, theta)
################################################################################


"""
    getindex(B::SSHBinner, theta, phi)

Returns the height of the bin which includes angles phi ∈ (0, 2pi) and theta ∈
(0, pi).
"""
function Base.getindex(B::AbstractSSH, theta::Real, phi::Real)
    0.0 <= phi <= 2pi || throw(ArgumentError(
        "Angle ϕ must be 0.0 ≤ ϕ = $phi ≤ π"
    ))
    0.0 <= theta <= pi || throw(ArgumentError(
        "Angle θ must be 0.0 ≤ θ = $theta ≤ π"
    ))

    # this is probably slow because the lambda function includes theta
    theta_index = findfirst(lower_bound -> theta + 1e-10 <= lower_bound, B.thetas)
    theta_index == nothing && (theta_index = length(B.thetas))
    theta_index == 1 || (theta_index -= 1)

    phi_index = trunc(Int64, phi * B.phi_N_over_2pi[theta_index]) + 1
    l = theta_index > 1 ? B.phi_cumsum[theta_index-1] : 0
    return get_value(B, l + phi_index)
end

# Colon phi
function Base.getindex(B::AbstractSSH, theta::Real, ::Colon)
    0.0 <= theta <= pi || throw(ArgumentError(
        "Angle θ must be 0.0 ≤ θ = $theta ≤ π"
    ))

    # this is probably slow because the lambda function includes theta
    theta_index = findfirst(lower_bound -> theta + 1e-10 <= lower_bound, B.thetas)
    theta_index == nothing && (theta_index = length(B.thetas))
    theta_index == 1 || (theta_index -= 1)

    l = theta_index > 1 ? B.phi_cumsum[theta_index-1] : 0
    return get_value.((B,), l .+ (1:B.phi_divisions[theta_index]))
end

# colon theta
function Base.getindex(B::AbstractSSH, ::Colon, phi::Real)
    0.0 <= phi <= 2pi || throw(ArgumentError(
        "Angle ϕ must be 0.0 ≤ ϕ = $phi ≤ π"
    ))

    idxs = Int64[]

    for theta_index in eachindex(B.phi_N_over_2pi)
        phi_index = trunc(Int64, phi * B.phi_N_over_2pi[theta_index]) + 1
        l = theta_index > 1 ? B.phi_cumsum[theta_index-1] : 0
        push!(idxs, l + phi_index)
    end

    return get_value.((B,), idxs)
end


# Array / Range + Colon/ Number
function Base.getindex(B::AbstractSSH, ::Colon, phis::AbstractArray)
    [getindex(B, phi, :) for phi in phis]
end
function Base.getindex(B::AbstractSSH, thetas::AbstractArray, ::Colon)
    [getindex(B, :, theta) for theta in thetas]
end
function Base.getindex(B::AbstractSSH, theta::Real, phis::AbstractArray)
    [getindex(B, phi, theta) for phi in phis]
end
function Base.getindex(B::AbstractSSH, thetas::AbstractArray, phi::Real)
    [getindex(B, phi, theta) for theta in thetas]
end
function Base.getindex(B::AbstractSSH, thetas::AbstractArray, phis::AbstractArray)
    [[getindex(B, phi, theta) for phi in phis] for theta in thetas]
end
function Base.getindex(B::AbstractSSH, ::Colon, ::Colon)
    last = 0
    out = Vector{Vector}(undef, length(B.phi_divisions))
    for (i, step) in enumerate(B.phi_divisions)
        out[i] = B.bins[last+1:last+step]
        last += step
    end
    out
end
