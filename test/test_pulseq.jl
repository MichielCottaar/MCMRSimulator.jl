@testset "test_pulseq.jl" begin
    directory = joinpath(pwd(), "example_pulseseq")
    @testset "read v1.3.1 fiddisp.seq file" begin
        seq = mr.read_pulseq(joinpath(directory, "fiddisp_v131.seq"))

        # starting with RF pulse
        @test length(seq.pulses) == 1
        pulse = seq.pulses[1]
        @test mr.start_time(pulse) == 0.1  # start after 100 microseconds delay
        @test mr.end_time(pulse) == 0.22  # pulse lasts 120 microseconds
        for cshape in (pulse.amplitude, pulse.phase)
            @test mr.start_time(cshape) ≈ 0.1
            @test mr.end_time(cshape) ≈ 0.22 
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
        @test length(seq.gradients) == 0

        # Single readout
        @test length(seq.readout_times) == 1
        rt = seq.readout_times[1]
        @test rt ≈ 0.22 + 5 + 0.02 + (1024 * 0.3125) / 2

        @test seq.TR ≈ 0.22 + 5 + 0.02 + (1024 * 0.3125)
    end
    @testset "read all v1.4.0 files in 01_from_FID_to_PRESS" begin
        path = joinpath(directory, "01_from_FID_to_PRESS_v140")
        for fn in readdir(path)
            if !endswith(fn, ".seq")
                continue
            end
            full_filename = joinpath(path, fn)
            seq = mr.read_pulseq(full_filename)
        end
    end
    @testset "check 01_from_FID_to_PRESS/01_FID.seq" begin
        fn = joinpath(directory, "01_from_FID_to_PRESS_v140", "01_FID.seq")
        seq = mr.read_pulseq(fn)
        @test length(seq.pulses) == 1
        @test mr.flip_angle(seq.pulses[1]) ≈ 90
        @test mr.start_time(seq.pulses[1]) ≈ 0.1  # determined by RF dead time
        @test mr.end_time(seq.pulses[1]) ≈ 0.6  # RF dead time + RF duration
        @test mr.phase(seq.pulses[1], 0.3) == 0
        @test length(seq.readout_times) == 1
        @test seq.readout_times[1] ≈ (
            0.6 +  # end of RF pulse
            0.02 +  # RF ringdown time (added by block duration)
            29.730 +  # Delay until ADC start (to get start at ADC at TE=30)
            256 / 2  # Halfway the ADC duration (256 ms)
        )
    end
    @testset "check 01_from_FID_to_PRESS/06_PRESS_center.seq" begin
        fn = joinpath(directory, "01_from_FID_to_PRESS_v140", "06_PRESS_center.seq")
        seq = mr.read_pulseq(fn)
        excitation = seq.pulses[1]
        @test length(seq.pulses) == 3
        #@test mr.flip_angle(seq.pulses[1]) ≈ 90
        @test mr.start_time(seq.pulses[1]) ≈ 0.1  # determined by RF dead time
        @test mr.end_time(seq.pulses[1]) ≈ 3.1  # RF dead time + RF duration
        @test mr.phase(seq.pulses[1], 0.3) ≈ 0 + rad2deg(0.5)
        @test mr.phase(seq.pulses[1], 1.6) ≈ 0
        @test length(seq.readout_times) == 1
        #@test seq.readout_times[1] ≈ 0.1 + 1.5 + 120 + 4096 * 62500 * 1e-6 / 2
        for p in seq.pulses[2:end]
            @test mr.end_time(p) - mr.start_time(p) ≈ 3  # should have been 4 for refocus pulses, but there is an error in the matlab generation
            @test mr.phase(p, mr.start_time(p)) ≈ 90 rtol=1e-5
            @test mr.phase(p, mr.start_time(p) + 0.38) ≈ 90 + rad2deg(0.5) rtol=1e-5
            @test mr.phase(p, mr.start_time(p) + 1.5) ≈ 90 rtol=1e-5
        end
    end
end