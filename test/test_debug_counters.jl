using Test
import MCMRSimulator.DebugCounters

abstract type DebugCounterProbe end
mutable struct GenericDebugCounterProbe <: DebugCounterProbe end
mutable struct SpecializedDebugCounterProbe <: DebugCounterProbe end
mutable struct UnregisteredDebugCounterProbe end

@testset "DebugCounters DebugTable" begin
    metrics = DebugCounters.DebugTable([
        (:cells_visited, "Cells visited per query", DebugCounters.FrequencyCounter(Int)),
        (:wide_integer, "Wide integer key", DebugCounters.FrequencyCounter(Int128)),
        (:net_count, "Net count", DebugCounters.Counter()),
        (:segment_length, "Segment length", DebugCounters.Histogram(0:1:10)),
        (:diffusion, "Diffusion summary", DebugCounters.SummaryMeasures()),
    ])

    @test DebugCounters.FrequencyTable === DebugCounters.FrequencyCounter
    @test_throws ArgumentError DebugCounters.snapshot(DebugCounters.FrequencyCounter(Int))
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

    DebugCounters.record!(metrics, :net_count, 5)
    DebugCounters.record!(metrics, :net_count, -2)
    counter_snapshot = DebugCounters.snapshot(metrics, :net_count)
    @test counter_snapshot isa DebugCounters.Counter{Int}
    @test counter_snapshot.count == 3
    counter_snapshot.count = 99
    @test DebugCounters.snapshot(metrics, :net_count).count == 3
    @test_throws ArgumentError DebugCounters.record!(metrics, :net_count, 1.5)

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
        (:counter, "counter", DebugCounters.Counter()),
        (:summary, "summary", DebugCounters.SummaryMeasures()),
        (:histogram, "histogram", DebugCounters.Histogram(0:1:10)),
    ])
    Threads.@threads for _ in 1:2000
        DebugCounters.record!(concurrent, :count, 1)
        DebugCounters.record!(concurrent, :counter, 2)
        DebugCounters.record!(concurrent, :summary, 1.0)
        DebugCounters.record!(concurrent, :histogram, 1.0)
    end
    concurrent_snapshot = DebugCounters.snapshot(concurrent)
    @test concurrent_snapshot[:count].counts == Dict(1 => 2000)
    @test concurrent_snapshot[:counter].count == 4000
    @test concurrent_snapshot[:summary].count == 2000
    @test concurrent_snapshot[:histogram].counts[2] == 2000

    output = sprint(DebugCounters.report, metrics)
    @test occursin("# cells_visited: Cells visited per query", output)
    @test occursin("underflow = 1", output)
    DebugCounters.reset!(metrics)
    @test isempty(DebugCounters.snapshot(metrics, :cells_visited).counts)
    @test DebugCounters.snapshot(metrics, :net_count).count == 0
    @test DebugCounters.snapshot(metrics, :segment_length).counts == zeros(Int, 10)
    @test DebugCounters.snapshot(metrics, :diffusion).count == 0

    @test_throws ArgumentError DebugCounters.DebugTable([
        (:duplicate, "first", DebugCounters.FrequencyCounter(Int)),
        (:duplicate, "second", DebugCounters.FrequencyCounter(Int)),
    ])
    @test_throws ArgumentError DebugCounters.Histogram([0, 1, 1])
    @test_throws ArgumentError DebugCounters.Histogram([0])
    overflow_counter = DebugCounters.Counter(Int8)
    overflow_counter.count = typemax(Int8)
    @test_throws OverflowError DebugCounters.record!(overflow_counter, 1)
end

@testset "DebugCounters instance registry" begin
    DebugCounters.reset!()
    DebugCounters.register!(DebugCounterProbe, [
        (:calls, "Calls", DebugCounters.FrequencyCounter(Int)),
        (:values, "Values", DebugCounters.SummaryMeasures()),
    ])
    DebugCounters.register!(SpecializedDebugCounterProbe, [
        (:specialized_calls, "Specialized calls", DebugCounters.FrequencyCounter(Int)),
    ])

    first_instance = GenericDebugCounterProbe()
    second_instance = GenericDebugCounterProbe()
    specialized_instance = SpecializedDebugCounterProbe()

    DebugCounters.record!(first_instance, :calls, 1)
    DebugCounters.record!(first_instance, :values, 2.5)
    DebugCounters.record!(second_instance, :calls, 2)
    DebugCounters.record!(specialized_instance, :specialized_calls, 7)

    @test DebugCounters.snapshot(first_instance)[:calls].counts == Dict(1 => 1)
    @test DebugCounters.snapshot(first_instance)[:values].count == 1
    @test DebugCounters.snapshot(second_instance)[:calls].counts == Dict(2 => 1)
    specialized_snapshot = DebugCounters.snapshot(specialized_instance)
    @test specialized_snapshot[:specialized_calls] isa DebugCounters.FrequencyCounter{Int}
    @test specialized_snapshot[:specialized_calls].counts == Dict(7 => 1)
    @test_throws KeyError DebugCounters.record!(specialized_instance, :calls, 1)
    @test_throws ArgumentError DebugCounters.snapshot(UnregisteredDebugCounterProbe())

    concurrent_instance = GenericDebugCounterProbe()
    Threads.@threads for _ in 1:1500
        DebugCounters.record!(concurrent_instance, :calls, 3)
    end
    @test DebugCounters.snapshot(concurrent_instance)[:calls].counts == Dict(3 => 1500)

    report_output = sprint(DebugCounters.report, first_instance)
    @test occursin("# calls: Calls", report_output)
    DebugCounters.reset!(first_instance)
    @test isempty(DebugCounters.snapshot(first_instance)[:calls].counts)

    # Global reset drops per-instance tables but leaves registrations available.
    DebugCounters.reset!()
    @test isempty(DebugCounters.snapshot(first_instance)[:calls].counts)
    DebugCounters.record!(first_instance, :calls, 9)
    @test DebugCounters.snapshot(first_instance)[:calls].counts == Dict(9 => 1)

    @test_throws ArgumentError DebugCounters.register!(DebugCounterProbe, [])
    DebugCounters.reset!()
end
