using BenchmarkTools
using SphereSurfaceHistogram
SSH = SphereSurfaceHistogram

# Needs to be defined, will be called on benchmark
const SUITE = BenchmarkGroup()



SUITE["Histogram creation"] = BenchmarkGroup()

suite = SUITE["Histogram creation"]["partition_sphere2"] = BenchmarkGroup()
for i in Int64[1e3, 1e4, 1e5, 1e6]
    suite["N = $i"] = @benchmarkable SSH.partition_sphere2($(4pi/i))
end


suite = SUITE["Histogram creation"]["partition_sphere_optim2"] = BenchmarkGroup()
for i in Int64[1e3, 1e4, 1e5, 1e6]
    suite["N = $i"] = @benchmarkable SSH.partition_sphere_optim2($(4pi/i))
end


suite = SUITE["Histogram creation"]["SSHBinner"] = BenchmarkGroup()
for i in Int64[1e3, 1e4, 1e5, 1e6]
    suite["N = $i"] = @benchmarkable SSHBinner($i)
end


suite = SUITE["Histogram creation"]["bin_many!"] = BenchmarkGroup()
for N_bins in Int64[1e2, 1e3, 1e4, 1e5]
    B = SSHBinner(N_bins)
    _suite = suite["N_bins = $N_bins"] = BenchmarkGroup()
    for N_spins in Int64[1e2, 1e3, 1e4, 1e5]
        spins = rand_spin(N_spins)

        _suite[
            "N_spins = $N_spins"
        ] = @benchmarkable bin_many!($B, $spins)
    end
end
