@testset "test_radio_frequency.jl: Using finite RF pulses" begin
    @testset "Constant RF pulse without off-resonance" begin
        for phase in (0, 30)
            rf = mr.constant_pulse(0, 9, 90, phase0=phase)
            @test rf.max_amplitude ≈ 1/36
            seq = mr.Sequence(pulses=[rf], TR=100.)
            sim = mr.Simulation(seq)
            signal = mr.signal(100, sim, 0:0.1:10)
            @test all(mr.propose_times(sim, 0, 9) .== 0:0.1:9)
            @test all(mr.propose_times(mr.Simulation(seq, rf_rotation=10), 0, 9) .== 0:1:9)
            increasing = signal[1:90]
            constant = signal[91:end]
            @test all(abs.(mr.longitudinal.(constant)) .<= 1e-8)
            @test all(mr.transverse.(constant) .≈ 100.)
            @test all(mr.longitudinal.(increasing) .≈ cosd.(0:89) .* 100.)
            @test all(mr.transverse.(increasing) .≈ sind.(0:89) .* 100.)
            @test all(mr.phase.(signal[2:end]) .≈ (phase + 90))
        end
    end
    @testset "Constant RF pulse with off-resonance" begin
        rf = mr.constant_pulse(0, 9, 90, off_resonance=1.)
        @test rf.max_amplitude ≈ sqrt(1 + (90 / (360 * 9))^2)
        seq = mr.Sequence(pulses=[rf], TR=100.)
        sim = mr.Simulation(seq, off_resonance=1)
        signal = mr.signal(100, sim, 0:0.1:10)
        increasing = signal[1:90]
        constant = signal[91:end]
        @test all(abs.(mr.longitudinal.(constant)) .<= 3)
        @test all(isapprox.(mr.transverse.(constant), 100., rtol=0.1))
        @test all(isapprox.(mr.longitudinal.(increasing), cosd.(0:89) .* 100., rtol=0.1))
        @test all(isapprox.(mr.transverse.(increasing), sind.(0:89) .* 100., rtol=0.1))
    end
end
