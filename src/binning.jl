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
    area(ϕ1, ϕ2, θ1, θ2) = (ϕ2 - ϕ1) * (cos(θ1) - cos(θ2))
    new_dA = dA
    N = round(Int64, 4pi/dA)
    for i in 1:max_iter
        thetas, dphis, phi_divs = partition_sphere2(new_dA)
        next_N = sum(phi_divs)
        if N == next_N
            @info "Finished after $i iterations"
            break
        else
            @info "$N -> $next_N"
            N = next_N
            new_dA = 4pi / N
        end
    end
    @info "Optimized dA = $new_dA ($N)"
    partition_sphere2(new_dA)
end


function bin_many!(bins, thetas, dphis, phi_divs, values)
    zs = cos.(thetas) .- -1e-20 # 1...-1
    skipped_over = cumsum(phi_divs)

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
        # yes it's the same!
        # @assert theta_index == binary_find(value[3], zs) "$value"

        phi_index = if phi_divs[theta_index] == 1
            1
        else
            # Normalize xy-component
            # y >= 0 for 0..pi, else y < 0
            # acos(x) well defined for 0..pi
            R = sqrt(value[1]^2 + value[2]^2)
            phi = if value[2] >= 0
                acos(value[1]/R)
            else
                2pi - acos(value[1]/R)
            end

            # min(⌊ϕ / Δϕ⌋ + 1, max_index)
            # TODO This is also taking a lot of time
            min(
                floor(Int64, phi / dphis[theta_index]) + 1,
                phi_divs[theta_index]
            )
        end

        # "Skipped" shells
        l = theta_index > 1 ? skipped_over[theta_index-1] : 0
        bins[l+phi_index] += 1
    end
    bins
end
