include("torun.jl")

using BenchmarkTools

suite = BenchmarkGroup()

for (key, sim) in pairs(simulations)
    suite[key] = @benchmarkable evolve(1000, $sim) evals=1
end

tune!(suite);