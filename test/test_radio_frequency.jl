@testset "Using finite RF pulses (test_radio_frequency.jl)" begin
    @testset "Constant RF pulse without off-resonance" begin
        for phase in (0, 30)
            rf = mr.constant_pulse(0, 9, 90, phase0=phase)
            @test rf.max_amplitude ≈ deg2rad(90) / 9
            seq = mr.Sequence(pulses=[rf], TR=100.)
            sim = mr.Simulation(seq)
            signal = mr.signal(100, sim, 0:0.1:10)
            @test all(mr.propose_times(sim, 0, 9) .== 0:0.1:9)
            increasing = signal[1:90]
            constant = signal[91:end]
            @test all(abs.(mr.longitudinal.(constant)) .<= 1e-8)
            @test all(mr.transverse.(constant) .≈ 100.)
            @test all(mr.longitudinal.(increasing) .≈ cos.(deg2rad.(0:89)) .* 100.)
            @test all(mr.transverse.(increasing) .≈ sin.(deg2rad.(0:89)) .* 100.)
            @test all(mr.phase.(signal[2:end]) .≈ (phase + 90))
        end
    end
    @testset "Constant RF pulse with off-resonance" begin
        rf = mr.constant_pulse(0, 9, 90, off_resonance=1.)
        @test rf.max_amplitude ≈ sqrt(1 + (deg2rad(90) / 9)^2)
        seq = mr.Sequence(pulses=[rf], TR=100.)
        sim = mr.Simulation(seq, off_resonance=1/(3. * mr.gyromagnetic_ratio))
        signal = mr.signal(100, sim, 0:0.1:10)
        increasing = signal[1:90]
        constant = signal[91:end]
        @test all(abs.(mr.longitudinal.(constant)) .<= 3)
        @test all(isapprox.(mr.transverse.(constant), 100., rtol=0.1))
        @test all(isapprox.(mr.longitudinal.(increasing), cos.(deg2rad.(0:89)) .* 100., rtol=0.1))
        @test all(isapprox.(mr.transverse.(increasing), sin.(deg2rad.(0:89)) .* 100., rtol=0.1))
    end
end
