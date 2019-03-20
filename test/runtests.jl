using SphereSurfaceHistogram
using Test

# Doesn't work here :(
# include(SphereSurfaceHistogram.PATHS[:point_gen])
# include(SphereSurfaceHistogram.PATHS[:Makie])
function to_cartesian(theta, phi)
    [
        sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)
    ]
end

function dual_points(B::SSHBinner)
    points = [[0., 0., 1.]]

    for i in 2:length(B.thetas)-2
        theta = 0.5(B.thetas[i] + B.thetas[i+1])
        N = B.phi_divisions[i]

        for k in 0:N-1
            phi = 2pi * (k + 0.5) / N
            push!(points, to_cartesian(theta, phi))
        end
    end
    push!(points, [0., 0., -1.])

    points
end

@testset "SphereSurfaceHistogram.jl" begin

    @testset "Flat Histogram" begin
        for _ in 1:10
            N = floor(Int64, 1_000_000rand()) + 1_000
            B = SSHBinner(N)
            bin_many!(B, dual_points(B))
            @test all(B.bins .== 1)
        end
    end
end
