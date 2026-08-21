include("torun.jl")

using BenchmarkTools

SUITE = BenchmarkGroup()

for (key, sim) in pairs(simulations)
    SUITE[key] = @benchmarkable readout(1000, $sim) evals=1
end

tune!(SUITE);
