using Test
import MRSimulator: Spin, Microstructure, evolve_to_time, time, field,
    gyromagnetic_ratio, RFPulse, apply_pulse, phase, longitudinal, transverse, time, position, 
    norm_angle, evolve, Sequence, relax
using StaticArrays

@testset "MRSimulator.jl" begin
    @testset "Generate and apply microstructural fields" begin
        @testset "Simulate empty environment/sequence" begin
            spin = evolve_to_time(Spin(), Microstructure(), 1.)
            @test time(spin) == 1.
            @test position(spin) == SA_F64[0., 0., 0.]
            @test phase(spin) == 0.
        end
        @testset "Test constant off-resonance field" begin
            spin = evolve_to_time(Spin(), Microstructure(off_resonance=field(2.)), 0.3, 3.)
            @test time(spin) == 0.3
            @test position(spin) == SA_F64[0., 0., 0.]
            @test phase(spin) ≈ norm_angle(rad2deg(0.6 * 3 * gyromagnetic_ratio))
        end
        @testset "Test gradient off-resonance field" begin
            micro = Microstructure(off_resonance=field(SA_F64[1.5, 0., 0.], 2.))
            spin = evolve_to_time(Spin(), micro, 0.3, 3.)
            @test time(spin) == 0.3
            @test position(spin) == SA_F64[0., 0., 0.]
            @test phase(spin) ≈ norm_angle(rad2deg(0.6 * 3 * gyromagnetic_ratio))
            # Move spin and evolve a bit further in time
            spin = evolve_to_time(Spin(0.3, SA_F64[3., 0., 0.], spin.orientation), micro, 0.5, 3.)
            @test time(spin) == 0.5
            @test position(spin) == SA_F64[3., 0., 0.]
            @test phase(spin) ≈ norm_angle(rad2deg((0.6 + (2. + 1.5 * 3.) * 0.2) * 3 * gyromagnetic_ratio) - 360)
        end
    end
    @testset "Simple relaxation" begin
        orient = Spin(transverse=1., longitudinal=0.).orientation
        pos = zero(SVector{3, Real})
        @testset "R2 relaxation" begin
            env = Microstructure(R2=field(2.))(pos)
            @test relax(orient, env, 0.3).transverse ≈ exp(-0.6)
        end
        @testset "R1 relaxation" begin
            env = Microstructure(R1=field(2.))(pos)
            @test relax(orient, env, 0.3).longitudinal ≈ 1 - exp(-0.6)
        end
    end
    @testset "Apply RF pulses" begin
        @testset "0 degree pulses should do nothing" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = RFPulse(0., 0., pulse_phase)
                for spin_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                    spin = Spin(phase=spin_phase, transverse=1.)
                    spin = apply_pulse(pulse, spin.orientation)
                    @test phase(spin) ≈ spin_phase
                    @test longitudinal(spin) ≈ 1.
                    @test transverse(spin) ≈ 1.
                end
            end
        end
        @testset "180 degree pulses should flip longitudinal" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = RFPulse(0., 180., pulse_phase)
                spin = Spin()
                @test longitudinal(spin) == 1.
                spin = apply_pulse(pulse, spin.orientation)
                @test longitudinal(spin) ≈ -1.
            end
        end
        @testset "90 degree pulses should eliminate longitudinal" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = RFPulse(0., 90., pulse_phase)
                spin = Spin()
                @test longitudinal(spin) == 1.
                spin = apply_pulse(pulse, spin.orientation)
                @test longitudinal(spin) ≈ 0. atol=1e-12
            end
        end
        @testset "Spins with same phase as pulse are unaffected by pulse" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                for flip_angle in (10, 90, 120, 180)
                    pulse = RFPulse(0., flip_angle, pulse_phase)

                    spin = Spin(longitudinal=0., transverse=1., phase=pulse_phase)
                    spin = apply_pulse(pulse, spin.orientation)
                    @test longitudinal(spin) ≈ 0. atol=1e-12
                    @test transverse(spin) ≈ 1.
                    @test phase(spin) ≈ pulse_phase
                end
            end
        end
        @testset "180 pulses flips phase around axis" begin
            for spin_phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 180, 0.)

                spin = Spin(longitudinal=0., transverse=1., phase=spin_phase)
                spin = apply_pulse(pulse, spin.orientation)
                @test longitudinal(spin) ≈ 0. atol=1e-12
                @test transverse(spin) ≈ 1.
                @test phase(spin) ≈ -spin_phase
            end
            for pulse_phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 180, pulse_phase)

                spin = Spin(longitudinal=0., transverse=1., phase=0.)
                spin = apply_pulse(pulse, spin.orientation)
                @test longitudinal(spin) ≈ 0. atol=1e-12
                @test transverse(spin) ≈ 1.
                @test phase(spin) ≈ 2 * pulse_phase
            end
        end
        @testset "90 pulses flips longitudinal spin into transverse plane" begin
            for pulse_phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 90, pulse_phase)
                spin = apply_pulse(pulse, Spin().orientation)
                @test longitudinal(spin) ≈ 0. atol=1e-12
                @test transverse(spin) ≈ 1.
                @test phase(spin) ≈ pulse_phase + 90
            end
        end
        @testset "90 pulses flips transverse spin into longitudinal plane" begin
            for pulse_phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 90, pulse_phase)
                spin_phase = (pulse_phase + 90)
                spin = Spin(longitudinal=0, transverse=1., phase=spin_phase)
                spin = apply_pulse(pulse, spin.orientation)
                @test longitudinal(spin) ≈ -1.
                @test transverse(spin) ≈ 0. atol=1e-12
            end
        end
    end
    @testset "Evolve a single spin fully" begin
        @testset "Empty environment and sequence" begin
            spins = evolve(Spin(), Microstructure(), Sequence(2.8), store_every=0.5)
            time = 0.
            for spin in spins
                @test spin.time == time
                time += 0.5
            end
            @test length(spins) == 6
            spins2 = evolve([Spin()], Microstructure(), Sequence(2.8), store_every=0.5)[:, 1]
            @test spins == spins2
        end
        @testset "Ensure data is stored at final TR" begin
            spins = evolve(Spin(), Microstructure(), Sequence(2.), store_every=0.5)
            @test length(spins) == 5
        end
    end
end
