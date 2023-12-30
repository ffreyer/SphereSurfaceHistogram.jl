struct SphereTessellationMap
    # positions of bins in theta
    # these may be irregular
    thetas::Vector{Float64}

    # number of divisions in phi direction. First bin is at phi = 0
    phi_divisions::Vector{Int64}

    # helpers for fast pushing
    zs::Vector{Float64}
    phi_N_over_2pi::Vector{Float64}
    phi_cumsum::Vector{Int64}
    N_thetas::Float64
    theta_factor::Float64
end

################################################################################
### Generic Constructor
################################################################################

function SphereTessellationMap(thetas::Vector{Float64}, phi_divisions::Vector{Int64})
    # Precalculate N/2pi for push!
    # use slightly more than 2pi to guarantee phi * phi_N_over_pi ∈ [0, 1)
    PI_PLUS = 2.0pi + 1e-14
    phi_N_over_2pi = phi_divisions ./ PI_PLUS

    # Also precalculate cumulative number of bins for each theta (from 0..pi)
    phi_cumsum = cumsum(phi_divisions)

    # Precalculate z(θ) for push!
    # move slightly to negative values to guarante value[3] >= minimum(zs)
    zs = cos.(thetas) #.- 1e-14 # 1...-1
    zs[end] -= 1e-14

    # also shift theta to be save with float precision issues
    thetas[1] -= 1e-14
    thetas[end] += 1e-14

    # for converting z values/angles to indices
    # to avoid conversions to Float64
    N_thetas = Float64(length(thetas)-1)
    theta_factor = Float64(length(thetas)-1) / (pi + 1e-14)

    return SphereTessellationMap(thetas, phi_divisions, zs, phi_N_over_2pi, phi_cumsum, N_thetas, theta_factor)
end


SphereTessellationMap(N_bins::Real; kwargs...) = SphereTessellationMap(round(Int64, N_bins); kwargs...)
function SphereTessellationMap(N_bins::Int64; method=partition_sphere2)
    thetas, phi_divisions = partition_sphere_optim(4pi/N_bins, method=method)
    return SphereTessellationMap(thetas, phi_divisions)
end

"""
    partition_sphere2(dA)

Partitions a sphere into bins of size `dA` where number of divisions in ϕ
(0..2π) is restricted to a power of 2.

Note that for a random `dA` this will not produce equally sized bins.
(Specifically bins at z ≈ 0 will be wrong.) The algorithm relies on
`partition_sphere_optim` to optimize `dA` to a value fitting on the sphere
surface which happens during the creation of an `SSHBinner`.
"""
function partition_sphere2(dA)
    # Use symmetry - only 0..pi/2 needed (± some delta if center of dA @ pi/2)

    # Starting at some point theta1 with width dphi,
    # computes theta2 for an area dA
    next_theta(dA, theta1, dphi) = acos(max(min(cos(theta1) - dA / dphi, 1.0), -1.0))

    # First area
    # This will now be our goal theta-width
    dtheta_goal = next_theta(dA, 0.0, 2pi)

    # Divide 2pi by phi_div, such that phi = 2pi/n with n integer
    phi_divs = Int64[1]
    thetas = Float64[0.0, dtheta_goal]

    # Currently radius of circle, should be sidelength of square
    # pi r^2 = x^2 => sqrt(pi)r = x
    dtheta_goal *= sqrt(pi)

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
    partition_sphere1(dA)

Partitions a sphere into bins of size `dA` where number of divisions in ϕ
(0..2π) is not restricted, i.e. it can be any integer.

Note that for a random `dA` this will not produce equally sized bins.
(Specifically bins at z ≈ 0 will be wrong.) The algorithm relies on
`partition_sphere_optim` to optimize `dA` to a value fitting on the sphere
surface which happens during the creation of an `SSHBinner`.
"""
function partition_sphere1(dA,)
    next_theta(dA, theta1, dphi) = acos(clamp(cos(theta1) - dA / dphi, -1.0, 1.0))

    # First area
    # This will now be our goal theta-width
    dtheta_goal = next_theta(dA, 0.0, 2pi)

    # Devide 2pi by phi_div, such that phi = 2pi/n with n integer
    phi_divs = Int64[1]
    thetas = Float64[0.0, dtheta_goal]

    # Currently radius of circle, should be sidelength of square
    # pi r^2 = x^2 => sqrt(pi)r = x
    dtheta_goal *= sqrt(pi)

    while true
        phi_div = last(phi_divs)
        theta1 = last(thetas)

        # Find best dtheta
        dtheta = next_theta(dA, theta1, 2pi / phi_div) - theta1
        best_dtheta = dtheta
        while true
            phi_div += 1
            dtheta = next_theta(dA, theta1, 2pi / phi_div) - theta1
            if abs(dtheta - dtheta_goal) < abs(best_dtheta - dtheta_goal)
                # New dtheta better? Keep it and keep trying
                best_dtheta = dtheta
            else
                # new dtheta worse? undo last phi_div increment and stop
                phi_div -= 1
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
    partition_sphere_optim(dA[; max_iter = 10, method = partition_sphere2])

Optimizes the sphere partitioning by modifying dA (~ the number of bins). After
a few iterations, dA becomes consistent for all bins.
"""
function partition_sphere_optim(dA; max_iter = 10, method=partition_sphere2)
    # NOTE this converges really fast

    new_dA = dA
    N = round(Int64, 4pi/dA)
    _iterations = 0
    for i in 1:max_iter
        _iterations = i
        thetas, phi_divs = method(new_dA)
        next_N = sum(phi_divs)
        if N == next_N
            break
        else
            N = next_N
            new_dA = 4pi / N
        end
    end

    method(new_dA)
end


################################################################################
### binning (push, append)
################################################################################


@fastmath @inline function fast_theta_index_approximation(x, l::Float64)
    # trunc faster than floor
    # Brackets cause llvm code to be `fmul` rather than `afoldl`, much
    # faster (~100x according to BenchmarkTools)
    # equivalent to
    # trunc(Int64, l*(0.5 - (7x + x^5 + x^9 + x^17 + x^31) / 22)
    # which is ≈ trunc(Int64, l*acos(x))
    trunc(Int64,
        l * (0.5 - x*(
            0.3181818181818182 + ((((x*x)*x)*x)*0.045454545454545456)*(
                1.0 + (((x*x)*x)*x) * (
                    1.0 + (((((((x*x)*x)*x)*x)*x)*x)*x)*(
                        1.0 + (((((((((((((x*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)*x)
                    )
                )
            )
        ))
    ) + 1
end

# Note:
# Tested for 10, 1000, 10_000_000 bins, each with 1_000_000 random spins pushed
@inline function fast_theta_index_search_z(tess::SphereTessellationMap, z)
    index = fast_theta_index_approximation(z, tess.N_thetas)
    if z > tess.zs[index]
        while true
            index -= 1
            (z <= tess.zs[index]) && (return index)
        end
    else
        while true
            index += 1
            (tess.zs[index] < z) && (return index-1)
        end
    end
    return index
end

@inline function fast_theta_index_search_theta(tess::SphereTessellationMap, theta)
    index = trunc(Int64, tess.theta_factor * theta) + 1
    if theta < tess.thetas[index]
        while true
            index -= 1
            (theta >= tess.thetas[index]) && (return index)
        end
    else
        while true
            index += 1
            (tess.thetas[index] > theta) && (return index-1)
        end
    end
    return index
end

function bin_index(tess::SphereTessellationMap, vec)
    theta_index = fast_theta_index_search_z(tess, vec[3])
    phi_index = if tess.phi_divisions[theta_index] == 1
        1
    else
        # Normalize xy-component
        # y >= 0 for 0..pi, else y < 0
        # acos(x) well defined for 0..pi
        R = sqrt(vec[1]^2 + vec[2]^2)
        phi = if vec[2] >= 0.0
            acos(vec[1]/R)
        else
            2.0pi - acos(vec[1]/R)
        end

        trunc(Int64, phi * tess.phi_N_over_2pi[theta_index]) + 1
    end

    # Calculate position in bins
    l = theta_index > 1 ? tess.phi_cumsum[theta_index-1] : 0
    return l + phi_index
end

function bin_index(tess::SphereTessellationMap, theta, phi)
    theta_index = fast_theta_index_search_theta(tess, theta)
    phi_index = trunc(Int64, tess.phi_N_over_2pi[theta_index] * phi) + 1

    # Calculate position in bins
    l = theta_index > 1 ? tess.phi_cumsum[theta_index-1] : 0
    return l + phi_index
end
