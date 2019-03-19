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

    return thetas, 2pi ./ phi_divs, phi_divs
end


# NOTE
# Wow this converges really fast
function partition_sphere_optim2(dA, max_iter = 10)
    # area(ϕ1, ϕ2, θ1, θ2) = (ϕ2 - ϕ1) * (cos(θ1) - cos(θ2))
    new_dA = dA
    N = round(Int64, 4pi/dA)
    _iterations = 0
    for i in 1:max_iter
        _iterations = i
        thetas, dphis, phi_divs = partition_sphere2(new_dA)
        next_N = sum(phi_divs)
        if N == next_N
            break
        else
            N = next_N
            new_dA = 4pi / N
        end
    end

    thetas, dphis, phi_divs = partition_sphere2(new_dA)

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

    thetas, dphis, phi_divs
end


const PI_PLUS = Float64(pi) + 1e-14
function bin_many!(bins, thetas, dphis, phi_divs, values)
    zs = cos.(thetas) .- 1e-10 # 1...-1
    skipped_over = cumsum(phi_divs)
    # TNOPI = (length(thetas) - 1) / PI_PLUS
    PNOPIs = phi_divs / 2PI_PLUS

    for value in values
        # 1...N-1
        # TODO Needs to be optimized - binary search tree?
        # welp...
        # theta_index = binary_find(value[3], zs)
        theta_index = 0
        for i in 1:length(zs)-1
            if value[3] >= zs[i+1]
                theta_index = i
                break
            end
        end

        # NOTE MUCH faster! ~28%!
        # NOTE but also totally wrong because thetas aren't linearly defined!
        # theta = acos(value[3])
        # theta_index = floor(Int64, theta * TNOPI) + 1
        # if theta_index != floor(Int64, theta * TNOPI) + 1
        #     println(
        #         "$theta_index != $(floor(Int64, theta * TNOPI) + 1)",
        #         " for $theta ($(value[3]))"
        #     )
        # end


        # yes it's the same!
        # @assert theta_index == binary_find(value[3], zs) "$value"

        phi_index = if phi_divs[theta_index] == 1
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

            # min(⌊ϕ / Δϕ⌋ + 1, max_index)
            # TODO This is also taking a lot of time
            # min(
            #     floor(Int64, phi / dphis[theta_index]) + 1,
            #     phi_divs[theta_index]
            # )
            # A bit faster ~ 5%
            floor(Int64, phi * PNOPIs[theta_index]) + 1
        end

        # "Skipped" shells
        l = theta_index > 1 ? skipped_over[theta_index-1] : 0
        bins[l+phi_index] += 1
    end
    bins
end
