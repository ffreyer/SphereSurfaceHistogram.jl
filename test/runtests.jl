using SphereSurfaceHistogram, GeometryBasics
using Test

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
        for method in (partition_sphere1, partition_sphere2)
            for _ in 1:10
                N = floor(Int64, 100_000rand()) + 1_000
                B = SSHBinner(N, method=method)

                # bin_positions (should) return center positions of bins in order
                # check if each point maps to the correct bin
                bin_centers = bin_positions(B)
                center_matched_bin = Bool[]
                for (i, v) in enumerate(bin_centers)
                    push!(B, v)
                    push!(center_matched_bin, B.bins[i] == 1)
                end
                @test all(center_matched_bin)

                # And again with angles
                empty!(B)
                empty!(center_matched_bin)
                idx = 1
                push!(B, 0.0, 0.0)
                push!(center_matched_bin, B.bins[idx] == 1)
                idx += 1
                for i in 2:length(B.tessellation.thetas)-2
                    theta = 0.5(B.tessellation.thetas[i] + B.tessellation.thetas[i+1])
                    N = B.tessellation.phi_divisions[i]

                    for k in 0:N-1
                        phi = 2pi * (k + 0.5) / N
                        push!(B, theta, phi)
                        push!(center_matched_bin, B.bins[idx] == 1)
                        idx += 1
                    end
                end
                push!(B, pi, 0.0)
                push!(center_matched_bin, B.bins[idx] == 1)
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
