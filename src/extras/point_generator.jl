
"""
    random_unit_vector(N)

Generates N random unit vectors, i.e. random points on a unit sphere.
"""
function random_unit_vector(N::Int64)
    phis = 2 * pi * rand(Float64, N)
    cts = 2 * rand(Float64, N) .- 1
    sts = sqrt.(1. .- cts .* cts) # sin(acos(cts)) # max(0., )

    [SVector{3, Float64}(
        sts[i] .* cos.(phis[i]),
        sts[i] .* sin.(phis[i]),
        cts[i]
    ) for i in 1:N]
end
