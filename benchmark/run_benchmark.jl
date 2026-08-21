using PkgBenchmark
using Markdown

mkconfig(; kwargs...) =
    BenchmarkConfig(
        env = Dict(
            "JULIA_NUM_THREADS" => "1",
            "OMP_NUM_THREADS" => "1",
        );
        kwargs...
    )

retune = get(ENV, "RETUNE", "false") == "true"

group_target = benchmarkpkg(
    dirname(@__DIR__),
    mkconfig(),
    retune=retune,
    resultfile = joinpath(@__DIR__, "result-target.json"),
)

export_markdown(stdout, group_target)

group_baseline = benchmarkpkg(
    dirname(@__DIR__),
    mkconfig(id = "add-benchmark"),
    retune=retune,
    resultfile = joinpath(@__DIR__, "result-baseline.json"),
)

judgement = judge(group_target, group_baseline)

export_markdown(stdout, judgement)
