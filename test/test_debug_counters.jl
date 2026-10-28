using Test
import MCMRSimulator.DebugCounters

@testset "DebugCounters DebugTable" begin
    metrics = DebugCounters.DebugTable([
        (:cells_visited, "Cells visited per query", DebugCounters.FrequencyCounter(Int)),
        (:wide_integer, "Wide integer key", DebugCounters.FrequencyCounter(Int128)),
        (:segment_length, "Segment length", DebugCounters.Histogram(0:1:10)),
        (:diffusion, "Diffusion summary", DebugCounters.SummaryMeasures()),
    ])

    @test DebugCounters.FrequencyTable === DebugCounters.FrequencyCounter
    @test_throws MethodError DebugCounters.snapshot(DebugCounters.FrequencyCounter(Int))
    DebugCounters.record!(metrics, :cells_visited, 0)
    DebugCounters.record!(metrics, :cells_visited, 2)
    foreach(_ -> DebugCounters.record!(metrics, :cells_visited, 2), 1:3)
    @test DebugCounters.snapshot(metrics, :cells_visited).counts == Dict(0 => 1, 2 => 4)

    wide_value = Int128(2)^100
    DebugCounters.record!(metrics, :wide_integer, wide_value)
    @test DebugCounters.snapshot(metrics, :wide_integer).counts == Dict(wide_value => 1)
    @test_throws ArgumentError DebugCounters.record!(metrics, :cells_visited, 1.5)
    @test_throws KeyError DebugCounters.record!(metrics, :missing, 1)
    @test_throws KeyError DebugCounters.snapshot(metrics, :missing)

    for value in (-1.0, 0.0, 0.5, 1.0, 9.9, 10.0, 11.0)
        DebugCounters.record!(metrics, :segment_length, value)
    end
    histogram = DebugCounters.snapshot(metrics, :segment_length)
    @test histogram.edges == collect(0.0:1.0:10.0)
    @test histogram.counts == [2, 1, 0, 0, 0, 0, 0, 0, 0, 2]
    @test histogram.underflow == 1
    @test histogram.overflow == 1

    for value in (1.0, 3.0, 5.0)
        DebugCounters.record!(metrics, :diffusion, value)
    end
    summary = DebugCounters.snapshot(metrics, :diffusion)
    @test summary.count == 3
    @test summary.minimum == 1.0
    @test summary.maximum == 5.0
    @test summary.mean == 3.0
    @test summary.m2 / (summary.count - 1) == 4.0

    # Table snapshots do not expose mutable collector storage.
    all_metrics = DebugCounters.snapshot(metrics)
    all_metrics[:cells_visited].counts[2] = 100
    all_metrics[:segment_length].counts[1] = 100
    @test DebugCounters.snapshot(metrics, :cells_visited).counts == Dict(0 => 1, 2 => 4)
    @test DebugCounters.snapshot(metrics, :segment_length).counts[1] == 2

    # A single table lock protects concurrent updates across its collectors.
    concurrent = DebugCounters.DebugTable([
        (:count, "count", DebugCounters.FrequencyCounter(Int)),
        (:summary, "summary", DebugCounters.SummaryMeasures()),
        (:histogram, "histogram", DebugCounters.Histogram(0:1:10)),
    ])
    Threads.@threads for _ in 1:2000
        DebugCounters.record!(concurrent, :count, 1)
        DebugCounters.record!(concurrent, :summary, 1.0)
        DebugCounters.record!(concurrent, :histogram, 1.0)
    end
    concurrent_snapshot = DebugCounters.snapshot(concurrent)
    @test concurrent_snapshot[:count].counts == Dict(1 => 2000)
    @test concurrent_snapshot[:summary].count == 2000
    @test concurrent_snapshot[:histogram].counts[2] == 2000

    output = sprint(DebugCounters.report, metrics)
    @test occursin("# cells_visited: Cells visited per query", output)
    @test occursin("underflow = 1", output)
    DebugCounters.reset!(metrics)
    @test isempty(DebugCounters.snapshot(metrics, :cells_visited).counts)
    @test DebugCounters.snapshot(metrics, :segment_length).counts == zeros(Int, 10)
    @test DebugCounters.snapshot(metrics, :diffusion).count == 0

    @test_throws ArgumentError DebugCounters.DebugTable([
        (:duplicate, "first", DebugCounters.FrequencyCounter(Int)),
        (:duplicate, "second", DebugCounters.FrequencyCounter(Int)),
    ])
    @test_throws ArgumentError DebugCounters.Histogram([0, 1, 1])
    @test_throws ArgumentError DebugCounters.Histogram([0])
end
