@testset "test_radio_frequency.jl: Using finite RF pulses" begin
    @testset "Constant RF pulse without off-resonance" begin
        for phase in (0, 30)
            rf = mr.constant_pulse(0, 9, 90, phase0=phase)
            @test rf.max_amplitude ≈ 1/36
            seq = mr.Sequence(components=[rf], TR=100.)
            sim = mr.Simulation(seq, rf_rotation=1.)
            signal = mr.readout(100, sim, 0:0.1:10)
            @test all(mr.propose_times(sim, 0, 9) .== 0:0.1:9)
            @test all(mr.propose_times(mr.Simulation(seq, rf_rotation=10), 0, 9) .== 0:1:9)
            increasing = signal[1:90]
            constant = signal[91:end]
            @test all(abs.(mr.longitudinal.(constant)) .<= 1e-8)
            @test all(mr.transverse.(constant) .≈ 100.)
            @test all(mr.longitudinal.(increasing) .≈ cosd.(0:89) .* 100.)
            @test all(mr.transverse.(increasing) .≈ sind.(0:89) .* 100.)
            @test all(mr.phase.(signal[2:end]) .≈ (phase - 90))
        end
    end
    @testset "Constant RF pulse with off-resonance" begin
        rf = mr.constant_pulse(0, 9, 90, off_resonance=1.)
        @test rf.max_amplitude ≈ sqrt(1 + (90 / (360 * 9))^2)
        seq = mr.Sequence(components=[rf], TR=100.)
        sim = mr.Simulation(seq, off_resonance=1, rf_rotation=1)
        signal = mr.readout(100, sim, 0:0.1:10)
        increasing = signal[1:90]
        constant = signal[91:end]
        @test all(abs.(mr.longitudinal.(constant)) .<= 3)
        @test all(isapprox.(mr.transverse.(constant), 100., rtol=0.1))
        @test all(isapprox.(mr.longitudinal.(increasing), cosd.(0:89) .* 100., rtol=0.1))
        @test all(isapprox.(mr.transverse.(increasing), sind.(0:89) .* 100., rtol=0.1))
    end
    @testset "Adiabatic pulse test" begin
        # Adopted from a jupyter notebook from Will Clarke
        ω = 0.52
        t = 31/2
        μ = 9.5
        f = 0.81
        tstep = 0.01

        t_axis = -t:tstep:t
    
        beta = f * π / μ
    
        amplitude = @. ω * (1 / cosh(beta * t_axis))
        phase = @. log(amplitude / ω) * μ

        pulse = mr.RFPulse(t_axis .+ t .+ 0.5, amplitude, rad2deg.(phase))
        grad = mr.MRGradients([0, 2 * t + 1], [1, 1])
        seq = mr.Sequence(components=[pulse, grad, mr.Readout(2 * t + 1)], TR=2 * t + 2)
        sim = mr.Simulation(seq, diffusivity=0.)

        get_signal(off_resonance) = mr.readout(mr.Spin(position=[off_resonance, 0, 0]), sim)
    
        @test mr.longitudinal(get_signal(0.)) <= -0.9
        @test mr.longitudinal(get_signal(1.)) >= 0.9
        @test mr.longitudinal(get_signal(-1.)) >= 0.9
        @test mr.transverse(get_signal(0.42)) >= 0.9
        @test mr.transverse(get_signal(-0.42)) >= 0.9
    end
    @testset "Excitation pulse" begin
        # Adopted from a jupyter notebook from Will Clarke
        t = 3
        f = 5.5
        b = 0.4
        tstep = 0.01
        t_axis = -t:tstep:t

        ampl = @. (
            π / 2 * 
            sin(π * t_axis * f) / (π * t_axis * f) *
            exp(-b^2 * t_axis^2)
            )
        ampl[t_axis .== 0] .= π / 2
        ampl .*= (90 / 102.79339283389233)

        pulse = mr.RFPulse(t_axis .+ (t + 1e-3), ampl)
        grad = mr.MRGradients([0, 2 * t + 1], [1, 1])
        seq = mr.Sequence(components=[pulse, grad, mr.Readout(2 * t + 1)], TR=2 * t + 2)
        sim = mr.Simulation(seq, diffusivity=0.)

        get_signal(off_resonance) = mr.readout(mr.Spin(position=[off_resonance, 0, 0]), sim)

        for o in (-2, 0, 2)
            @test abs(mr.longitudinal(get_signal(o))) <= 0.3
            @test mr.transverse(get_signal(o)) >= 0.9
        end

        for o in (-10, -5, 5, 10)
            @test mr.longitudinal(get_signal(o)) >= 0.9
            @test mr.transverse(get_signal(o)) <= 0.2
        end

        @testset "Spins with very short T2" begin
            sim1 = mr.Simulation(seq, diffusivity=0., R2=100, rf_rotation=5) # T2 of 10 ns
            sim2 = mr.Simulation(seq, diffusivity=0., R2=100, rf_rotation=0.5)

            for o in (0, 5, 20, 50)
                spin = mr.Spin(position=[o, 0, 0])
                s1 = mr.readout(spin, sim1)
                s2 = mr.readout(spin, sim2)
                @test (1 - mr.longitudinal(s1)) ≈ (1 - mr.longitudinal(s2)) rtol=0.01
                @test mr.transverse(s1) < 1e-3
                @test mr.transverse(s2) < 1e-3
            end
        end
    end
end
