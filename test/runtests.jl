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

    # These should pass without issues
    # If they don't the binning fails for these edge cases.
    @testset "Edge Cases" begin
        B = SSHBinner(1_000)
        push!(B, [0., 0., 1.])
        @test true
        push!(B, [0., 0., -1.])
        @test true
        push!(B, [0., 1., 0.])
        @test true
        push!(B, [0., -1., 0.])
        @test true
        push!(B, [1., 0., 0.])
        @test true
        push!(B, [-1., 0., 0.])
        @test true
    end


    @testset "Flat Histogram" begin
        for method in [
            SphereSurfaceHistogram.partition_sphere1,
            SphereSurfaceHistogram.partition_sphere2
        ]
            for _ in 1:10
                N = floor(Int64, 1_000_000rand()) + 1_000
                B = SSHBinner(N, method=method)

                # dual_points (should) return center positions of bins in order
                # check if each point maps to the correct bin
                bin_centers = dual_points(B)
                center_matched_bin = Bool[]
                for (i, v) in enumerate(bin_centers)
                    push!(B, v)
                    push!(center_matched_bin, B.bins[i] == 1)
                end
                @test all(center_matched_bin)
            end
        end
    end

    @testset "Indexing" begin
        B = SSHBinner(10_000)
        append!(B, [[cos(pi), sin(pi), 0.0] for _ in 1:1000])

        # true position is 4969
        @test B.bins[4969] == 1000

        # indexing
        @test B[pi/2, pi] == 1000
        @test sum(B[pi/2, :]) == 1000
        @test sum(B[:, pi]) == 1000
    end
end
