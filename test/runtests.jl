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
        for T in (SSHBinner, SSHAverager)
            B = T(10_000)
            if T == SSHBinner
                append!(B, [[cos(pi), sin(pi), 0.0] for _ in 1:1000])
                # true position is 4969
                @test B.bins[4969] == 1000
            else
                append!(B, [[cos(pi), sin(pi), 0.0] for _ in 1:1000], ones(1000))
                # true position is 4969
                @test B.counts[4969] == 1000
                @test B.sums[4969] == 1000.0
            end

            result = SphereSurfaceHistogram.get_value(B, 4969)
            # indexing
            @test B[pi/2, pi] == result
            @test sum(B[pi/2, :]) == result
            @test sum(B[:, pi]) == result
        end
    end

    @testset "SSHAverager + SSHBinner utils" begin
        for T in (SSHBinner, SSHAverager)
            B = T(100, method = partition_sphere1)
            @test length(B) == 108
            if T == SSHBinner
                push!(B, 0.5, 0.2)
            else
                push!(B, 0.5, 0.2, 1.0)
            end
            @test minimum(B) ≈ 0.0
            @test maximum(B) ≈ 1.0
            # may be changed...
            @test size(B) == [1, 7, 12, 16, 18, 18, 16, 12, 7, 1]

            if T == SSHBinner
                append!(B, rand(100), rand(100))
                @test get_values(B) == B.bins
                C = T(B.bins, method = partition_sphere1)
                @test C == B
            else
                append!(B, rand(100), rand(100), rand(100))
                @test get_values(B) == B.sums ./ max.(1, B.counts)
                C = T(B.sums, B.counts, method = partition_sphere1)
                @test C == B
            end
        end

        for T in (SSHBinner, SSHAverager)
            B = T(100, method = partition_sphere2)
            @test length(B) == 98
            if T == SSHBinner
                push!(B, 0.5, 0.2)
            else
                push!(B, 0.5, 0.2, 1.0)
            end
            @test minimum(B) ≈ 0.0
            @test maximum(B) ≈ 1.0
            # may be changed...
            @test size(B) == [1, 8, 16, 16, 16, 16, 16, 8, 1]

            if T == SSHBinner
                append!(B, rand(100), rand(100))
                @test get_values(B) == B.bins
                C = T(B.bins, method = partition_sphere2)
                @test C == B
            else
                append!(B, rand(100), rand(100), rand(100))
                @test get_values(B) == B.sums ./ max.(1, B.counts)
                C = T(B.sums, B.counts, method = partition_sphere2)
                @test C == B
            end
        end
    end
end
