@testset "test_pulseseq.jl" begin
    directory = joinpath(pwd(), "example_pulseseq")
    @testset "read fiddisp.seq" begin
        keywords = open(mr.read_pulseseq_sections, joinpath(directory, "fiddisp.seq"))
        seq = mr.build_sequence(;keywords...);

        # starting with RF pulse
        @test length(seq.pulses) == 1
        pulse = seq.pulses[1]
        @show pulse
        @test mr.start_time(pulse) == 0.1  # start after 100 microseconds delay
        @test mr.end_time(pulse) == 0.22  # pulse lasts 120 microseconds
        for cshape in (pulse.amplitude, pulse.phase)
            @test mr.start_time(cshape) == 0.1
            @test mr.end_time(cshape) == 0.22 
        end
        for (time, ampl) in [
            (0.09, 0.),
            (0.11, 2.5),
            (0.19, 2.5),
            (0.21, 0.),
            (0.23, 0.),
        ]
            @test iszero(mr.phase(pulse, time))
            @test mr.amplitude(pulse, time) ≈ ampl
        end

        # no gradients
        @test mr.control_points(seq.gradient) == [0, seq.TR]

        # Single readout
        @test length(seq.readout_times) == 1
        rt = seq.readout_times[1]
        @assert rt ≈ 0.22 + 5 + 0.02 + (1024 * 0.3125) / 2

        @assert seq.TR ≈ 0.22 + 5 + 0.02 + (1024 * 0.3125)
    end
end