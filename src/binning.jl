struct SSHBinner
    thetas::Vector{Float64}
    phi_divisions::Vector{Int64}

    # bin_many helpers
    zs::Vector{Float64}
    phi_N_over_2pi::Vector{Float64}
    cumsum::Vector{Int64}

    bins::Vector{Int64}
end

function SSHBinner(thetas::Vector{Float64}, phi_divisions::Vector{Int64})
    # Precalculate N/2pi for bin_many!
    # use slightly more than 2pi to guarantee phi * phi_N_over_pi ∈ [0, 1)
    PI_PLUS = 2.0pi + 1e-14
    phi_N_over_2pi = phi_divisions ./ PI_PLUS

    # Precalculate z(θ) for bin_many!
    # move slightly to negative values to guarantee that for an input value
    # value[3] >= minimum(zs)
    zs = cos.(thetas) .- 1e-14 # 1...-1

    # Also precalculate cumulative number of bins for each theta (from 0..pi)
    _cumsum = cumsum(phi_divisions)

    N_bins = sum(phi_divisions)
    bins = zeros(Int64, N_bins)
    SSHBinner(thetas, phi_divisions, zs, phi_N_over_2pi, _cumsum, bins)
end

"""
    SSHBinner(N_bins)

Generate a Sphere Surface Histogram Binner with approximately N_bins.
"""
SSHBinner(N_bins::Real) = SSHBinner(round(Int64, N_bins))
function SSHBinner(N_bins::Int64)
    thetas, phi_divisions = partition_sphere_optim2(4pi/N_bins)
    SSHBinner(thetas, phi_divisions)
end


"""
    partition_sphere2(dA[, factor=sqrt(pi)])

Partitions a sphere into bins, where each bin except the ones around z≈0 has
equal area dA. The number of divisions in ϕ (0..2π) is restricted to powers of
2 and the divisions of θ vary to keep dA constant.
"""
function partition_sphere2(dA, factor=sqrt(pi))
    # Use symmetry - only 0..pi/2 needed (± some delta if center of dA @ pi/2)

    # Starting at some point theta1 with width dphi,
    # computes theta2 for an area dA
    next_theta(dA, theta1, dphi) = acos(max(min(cos(theta1) - dA / dphi, 1.0), -1.0))

    # First area
    # This will now be our goal theta-width
    dtheta_goal = next_theta(dA, 0.0, 2pi)

    # Devide 2pi by phi_div, such that phi = 2pi/n with n integer
    phi_divs = Int64[1]
    thetas = Float64[0.0, dtheta_goal]

    # Currently radius of circle, should be sidelength of square
    # pi r^2 = x^2 => sqrt(pi)r = x
    dtheta_goal *= factor

    while true
        phi_div = last(phi_divs)
        theta1 = last(thetas)

        # Find best dtheta
        dtheta = next_theta(dA, theta1, 2pi / phi_div) - theta1
        best_dtheta = dtheta
        while true
            phi_div *= 2
            dtheta = next_theta(dA, theta1, 2pi / phi_div) - theta1
            if abs(dtheta - dtheta_goal) < abs(best_dtheta - dtheta_goal)
                # New dtheta better? Keep it and keep trying
                best_dtheta = dtheta
            else
                # new dtheta worse? undo last phi_div increment and stop
                phi_div = div(phi_div, 2)
                break
            end
        end

        if (theta1 + best_dtheta / 2.0) > pi/2
            break
        else
            push!(thetas, theta1 + best_dtheta)
            push!(phi_divs, phi_div)
        end
    end

    dtheta = thetas[end] - thetas[end-1]
    theta_fourth = thetas[end-1] + 3dtheta/4
    if theta_fourth > pi/2
        # center should be flush
        thetas[end] = pi - thetas[end-1]
        for i in length(thetas)-2:-1:1
            push!(thetas, pi - thetas[i])
            push!(phi_divs, phi_divs[i])
        end
    else
        # last(thetas) should be flush
        thetas[end] = pi/2
        for i in length(thetas)-1:-1:1
            push!(thetas, pi - thetas[i])
            push!(phi_divs, phi_divs[i])
        end
    end

    return thetas, phi_divs
end


"""
    partition_sphere_optim2(dA[max_iter = 10])

Optimizes the sphere partitioning by modifying dA (the number of bins). After a
few iterations, dA becomes consistent for all bins.
"""
function partition_sphere_optim2(dA, max_iter = 10)
    # NOTE this converges really fast
    # area(ϕ1, ϕ2, θ1, θ2) = (ϕ2 - ϕ1) * (cos(θ1) - cos(θ2))
    new_dA = dA
    N = round(Int64, 4pi/dA)
    _iterations = 0
    for i in 1:max_iter
        _iterations = i
        thetas, phi_divs = partition_sphere2(new_dA)
        next_N = sum(phi_divs)
        if N == next_N
            break
        else
            N = next_N
            new_dA = 4pi / N
        end
    end

    thetas, phi_divs = partition_sphere2(new_dA)

    # Some stats
    # dAs = map(enumerate(dphis)) do t
    #     i, dphi = t
    #     dphi * (cos(thetas[i]) - cos(thetas[i+1]))
    # end
    # dA_mean = round(mean(dAs), sigdigits = 3)
    # dA_std = round(std(dAs), sigdigits = 3)
    # dA_perc = round(100dA_std/dA_mean, sigdigits = 3)
    # @info(
    #     "Optimized dA = $dA_mean ± $dA_std  ($(dA_perc)% @ " *
    #     "$(sum(phi_divs)) bins after $_iterations Iterations)"
    # )

    thetas, phi_divs
end


"""
    push!(B::SSHBinner, value)

Bins a single value (three dimensional unit vector).
"""
function Base.push!(B::SSHBinner, value)
    theta_index = 0
    for i in 1:length(B.zs)-1
        if value[3] >= B.zs[i+1]
            theta_index = i
            break
        end
    end

    phi_index = if B.phi_divisions[theta_index] == 1
        1
    else
        # Normalize xy-component
        # y >= 0 for 0..pi, else y < 0
        # acos(x) well defined for 0..pi
        R = sqrt(value[1]^2 + value[2]^2)
        phi = if value[2] >= 0.0
            acos(value[1]/R)
        else
            2pi - acos(value[1]/R)
        end

        floor(Int64, phi * B.phi_N_over_2pi[theta_index]) + 1
    end

    # Calculate position in bins
    l = theta_index > 1 ? B.cumsum[theta_index-1] : 0
    B.bins[l+phi_index] += 1
    nothing
end


"""
    append!(B::SSHBinner, values)

Bins an array of values (three dimensional unit vectors).
"""
function Base.append!(B::SSHBinner, values)
    for v in values
        push!(B, v)
    end
    nothing
end


################################################################################
# Cosmetics

function Base.show(io::IO, B::SSHBinner)
    print(io, "SSHBinner with $(length(B.bins)) bins")
end
