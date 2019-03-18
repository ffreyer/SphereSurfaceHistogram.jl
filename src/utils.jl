function binary_find(x, xs, start=1, stop=length(xs))
    width = stop - start + 1
    if width < 10
        theta_index = 0
        for i in start:stop-1
            if x >= xs[i+1]
                return i
            end
        end
    else
        middle = start + div(width, 2)
        if x <= xs[middle]
            return binary_find(x, xs, middle, stop)
        else
            return binary_find(x, xs, start, middle)
        end
    end
end


function rand_spin(N::Int64)
    phis = 2 * pi * rand(Float64, N)
    cts = 2 * rand(Float64, N) .- 1
    sts = sqrt.(1. .- cts .* cts) # sin(acos(cts)) # max(0., )

    [Point3{Float64}(
        sts[i] .* cos.(phis[i]),
        sts[i] .* sin.(phis[i]),
        cts[i]
    ) for i in 1:N]
end


function to_cartesian(theta, phi)
    Vec3f0(
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    )
end
