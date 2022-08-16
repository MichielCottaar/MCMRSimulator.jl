using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using PkgBenchmark
branch_or_commit = nothing
kwargs = Dict(
    :juliacmd => `julia -O3 --project=MRSimulator`, 
    :env => Dict("JULIA_NUM_THREADS" => 4))

if length(ARGS) > 0
    branch_or_commit = ARGS[1]
    result = judge(
        "MRSimulator",
        BenchmarkConfig(;kwargs...),
        BenchmarkConfig(;id=branch_or_commit, kwargs...)
    )
else
    result = benchmarkpkg(
        "MRSimulator",
        BenchmarkConfig(;kwargs...),
    )
end

export_markdown(stdout, result)