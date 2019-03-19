using BenchmarkTools
using SphereSurfaceHistogram
SSH = SphereSurfaceHistogram

# Needs to be defined, will be called on benchmark
const SUITE = BenchmarkGroup()



SUITE["Histogram creation"] = BenchmarkGroup()

suite = SUITE["Histogram creation"]["partition_sphere2"] = BenchmarkGroup()
for i in Int64[1e3, 1e4, 1e5, 1e6]
    suite["N = $i"] = @benchmarkable partition_sphere2($(4pi/i))
end


suite = SUITE["Histogram creation"]["partition_sphere_optim2"] = BenchmarkGroup()
for i in Int64[1e3, 1e4, 1e5, 1e6]
    suite["N = $i"] = @benchmarkable partition_sphere_optim2($(4pi/i))
end


suite = SUITE["Histogram creation"]["bin_many!"] = BenchmarkGroup()
for N_bins in Int64[1e3, 1e4, 1e5, 1e6]
    thetas, dphis, phi_divs = partition_sphere_optim2(4pi/N_bins)
    N = sum(phi_divs)
    bins = zeros(Int64, N)
    for N_spins in Int64[1e3, 1e4, 1e5, 1e6]
        spins = rand_spin(N_spins)

        suite[
            "N_bins = $N_bins; N_spins = $N_spins"
        ] = @benchmarkable bin_many!($bins, $thetas, $dphis, $phi_divs, $spins)
    end
end
