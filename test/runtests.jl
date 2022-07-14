using Test
import MRSimulator: Spin, ZeroField, Microstructure, evolve!, time, ConstantField, GradientField, gyromagnetic_ratio, RFPulse, apply!
using StaticArrays

@testset "MRSimulator.jl" begin
    @testset "Generate and apply microstructural fields" begin
        @testset "Simulate empty environment/sequence" begin
            spin = evolve!(Spin(), Microstructure(), 1.)
            @test spin.time == 1.
            @test spin.position == SA_F64[0., 0., 0.]
            @test spin.phase == 0.
        end
        @testset "Test constant off-resonance field" begin
            spin = evolve!(Spin(), Microstructure(off_resonance=ConstantField(2.)), 0.3, 3.)
            @test spin.time == 0.3
            @test spin.position == SA_F64[0., 0., 0.]
            @test spin.phase ≈ 0.6 * 3 * gyromagnetic_ratio
        end
        @testset "Test gradient off-resonance field" begin
            micro = Microstructure(off_resonance=GradientField(SA_F64[1.5, 0., 0.], 2.))
            spin = Spin()
            evolve!(spin, micro, 0.3, 3.)
            @test spin.time == 0.3
            @test spin.position == SA_F64[0., 0., 0.]
            @test spin.phase ≈ 0.6 * 3 * gyromagnetic_ratio
            # Move spin and evolve a bit further in time
            spin.position = SA_F64[3., 0., 0.]
            evolve!(spin, micro, 0.5, 3.)
            @test spin.time == 0.5
            @test spin.position == SA_F64[3., 0., 0.]
            @test spin.phase ≈ (0.6 + (2. + 1.5 * 3.) * 0.2) * 3 * gyromagnetic_ratio
        end
    end
    @testset "Apply RF pulses" begin
        @testset "0 degree pulses should do nothing" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = RFPulse(0., 0., pulse_phase)
                for spin_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                    spin = Spin(phase=spin_phase / 180 * π, transverse=1.)
                    apply!(pulse, spin)
                    @test spin.phase ≈ (spin_phase % 360) / 180 * π
                    @test spin.longitudinal ≈ 1.
                    @test spin.transverse ≈ 1.
                end
            end
        end
        @testset "180 degree pulses should flip longitudinal" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = RFPulse(0., 180., pulse_phase)
                spin = Spin()
                @test spin.longitudinal == 1.
                apply!(pulse, spin)
                @test spin.longitudinal ≈ -1.
            end
        end
        @testset "90 degree pulses should eliminate longitudinal" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = RFPulse(0., 90., pulse_phase)
                spin = Spin()
                @test spin.longitudinal == 1.
                apply!(pulse, spin)
                @test spin.longitudinal ≈ 0. atol=1e-12
            end
        end
        @testset "Spins with same phase as pulse are unaffected by pulse" begin
            for phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                for flip_angle in (10, 90, 120, 180)
                    pulse = RFPulse(0., flip_angle, phase)

                    spin_phase = phase / 180 * π
                    spin = Spin(longitudinal=0., transverse=1., phase=spin_phase)
                    apply!(pulse, spin)
                    @test spin.longitudinal ≈ 0. atol=1e-12
                    @test spin.transverse ≈ 1.
                    @test spin.phase ≈ spin_phase
                end
            end
        end
        @testset "180 pulses flips phase around axis" begin
            for phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 180, 0.)

                spin_phase = phase / 180 * π
                spin = Spin(longitudinal=0., transverse=1., phase=spin_phase)
                apply!(pulse, spin)
                @test spin.longitudinal ≈ 0. atol=1e-12
                @test spin.transverse ≈ 1.
                @test spin.phase ≈ -spin_phase
            end
            for phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 180, phase)

                spin_phase = phase / 180 * π
                spin = Spin(longitudinal=0., transverse=1., phase=0.)
                apply!(pulse, spin)
                @test spin.longitudinal ≈ 0. atol=1e-12
                @test spin.transverse ≈ 1.
                @test spin.phase ≈ 2 * phase / 180 * π
            end
        end
        @testset "90 pulses flips longitudinal spin into transverse plane" begin
            for phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 90, phase)
                spin = Spin()
                apply!(pulse, spin)
                @test spin.longitudinal ≈ 0. atol=1e-12
                @test spin.transverse ≈ 1.
                @test spin.phase ≈ (phase + 90) / 180 * π
            end
        end
        @testset "90 pulses flips transverse spin into longitudinal plane" begin
            for phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 90, phase)
                spin_phase = (phase + 90) / 180 * π
                spin = Spin(longitudinal=0, transverse=1., phase=spin_phase)
                apply!(pulse, spin)
                @test spin.longitudinal ≈ -1.
                @test spin.transverse ≈ 0. atol=1e-12
            end
        end
    end
end
