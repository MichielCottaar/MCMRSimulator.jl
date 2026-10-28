using Test
import MCMRSimulator.DebugCounters

@testset "DebugCounters" begin
    DebugCounters.reset!()
    @test DebugCounters.register!(:cells_visited, Int) isa DebugCounters.FrequencyTable
    DebugCounters.record!(:cells_visited, 0)
    DebugCounters.record!(:cells_visited, 2)
    foreach(_ -> DebugCounters.record!(:cells_visited, 2), 1:3)
    @test DebugCounters.snapshot()[:cells_visited] == Dict(0 => 1, 2 => 4)
    @test DebugCounters.register!(:cells_visited, Int) isa DebugCounters.FrequencyTable
    @test_throws ArgumentError DebugCounters.register!(:cells_visited, Float64)
    @test_throws ArgumentError DebugCounters.record!(:unknown, 1)
    @test_throws ArgumentError DebugCounters.record!(:cells_visited, 1.5)
    DebugCounters.register!(:wide_integer, Int128)
    wide_value = Int128(2)^100
    DebugCounters.record!(:wide_integer, wide_value)
    @test DebugCounters.snapshot()[:wide_integer] == Dict(wide_value => 1)

    histogram = DebugCounters.register!(:segment_length, Float64; histogram=true, bins=0:1:10)
    @test histogram isa DebugCounters.Histogram
    for value in (-1.0, 0.0, 0.5, 1.0, 9.9, 10.0, 11.0)
        DebugCounters.record!(:segment_length, value)
    end
    histogram_snapshot = DebugCounters.snapshot()[:segment_length]
    @test histogram_snapshot.edges == collect(0.0:1.0:10.0)
    @test histogram_snapshot.counts == [2, 1, 0, 0, 0, 0, 0, 0, 0, 2]
    @test histogram_snapshot.underflow == 1
    @test histogram_snapshot.overflow == 1

    # Snapshots are independent, and reset clears registrations and values.
    histogram_snapshot.counts[1] = 100
    @test DebugCounters.snapshot()[:segment_length].counts[1] == 2
    DebugCounters.reset!()
    @test isempty(DebugCounters.snapshot())

    @test_throws ArgumentError DebugCounters.register!(:bad_histogram, Float64; histogram=true)
    @test_throws ArgumentError DebugCounters.register!(:bad_edges, Int; histogram=true, bins=[0, 1, 1])
    @test_throws ArgumentError DebugCounters.register!(:unexpected_bins, Int; bins=0:1:10)
    @test_throws ArgumentError DebugCounters.register!(:abstract_type, Real)
    DebugCounters.register!(:float32, Float32)
    @test DebugCounters.register!(:float32, Float32) isa DebugCounters.SummaryMeasures
    DebugCounters.record!(:float32, Float32(0.5))
    @test_throws ArgumentError DebugCounters.record!(:float32, 0.1)

    DebugCounters.register!(:candidate_tests, Int)
    DebugCounters.record!(:candidate_tests, 1)
    output = sprint(DebugCounters.report)
    @test occursin("candidate_tests", output)
    @test occursin("1 => 1", output)
    DebugCounters.reset!()
end
