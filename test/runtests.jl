using Test
import MRSimulator: MRSimulator, Spin, Microstructure, evolve_to_time, time, field,
    gyromagnetic_ratio, RFPulse, apply_pulse, phase, longitudinal, transverse, time, position, 
    norm_angle, evolve, Sequence, relax, vector2spin, vector
using StaticArrays

@testset "MRSimulator.jl" begin
    @testset "Generate and apply microstructural fields" begin
        @testset "Simulate empty environment/sequence" begin
            spin = evolve_to_time(Spin(), Microstructure(), 0., 1.)
            @test position(spin) == SA_F64[0., 0., 0.]
            @test phase(spin) == 0.
        end
        @testset "Test constant off-resonance field" begin
            spin = evolve_to_time(Spin(), Microstructure(off_resonance=field(2.)), 0., 0.3, 3.)
            @test position(spin) == SA_F64[0., 0., 0.]
            @test phase(spin) ≈ norm_angle(rad2deg(0.6 * 3 * gyromagnetic_ratio))
        end
        @testset "Test gradient off-resonance field" begin
            micro = Microstructure(off_resonance=field(SA_F64[1.5, 0., 0.], 2.))
            spin = evolve_to_time(Spin(), micro, 0., 0.3)
            @test position(spin) == SA_F64[0., 0., 0.]
            @test phase(spin) ≈ norm_angle(rad2deg(0.6 * 3 * gyromagnetic_ratio))
            # Move spin and evolve a bit further in time
            spin = evolve_to_time(Spin(SA_F64[3., 0., 0.], spin.orientation), micro, 0.3, 0.5, 3.)
            @test position(spin) == SA_F64[3., 0., 0.]
            @test phase(spin) ≈ norm_angle(rad2deg((0.6 + (2. + 1.5 * 3.) * 0.2) * 3 * gyromagnetic_ratio))
        end
        @testset "Fields with different types" begin
            pos = SA_F64[1., 0., 0.]
            @test isa(field()(pos), Float64)
            @test isa(field(Int)(pos), Int)
            @test isa(field(0)(pos), Int)
            @test isa(field(0), MRSimulator.ZeroField{Int})
            @test isa(field(2)(pos), Int)
            @test isa(field(2), MRSimulator.ConstantField{Int})
            @test isa(field([0, 0, 0], 0), MRSimulator.ZeroField{Int})
            @test isa(field([0, 0, 0], 0)(pos), Int)
            @test isa(field([0, 0, 0], 2), MRSimulator.ConstantField{Int})
            @test isa(field([0, 0, 0], 2)(pos), Int)
            @test isa(field([1, 0, 0], 2), MRSimulator.GradientField{Int})
            @test isa(field([1, 0, 0], 2)(pos), Float64)
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
    @testset "Spin conversions" begin
        vec = SA_F64[1, 0, 0]
        @test vector2spin(vec).longitudinal ≈ 0.
        @test vector2spin(vec).transverse ≈ 1.
        @test phase(vector2spin(vec)) ≈ 0.
        @test vector(vector2spin(vec)) ≈ vec

        vec = SA_F64[0, 1, 1]
        @test vector2spin(vec).longitudinal ≈ 1.
        @test vector2spin(vec).transverse ≈ 1.
        @test phase(vector2spin(vec)) ≈ 90.
        @test vector(vector2spin(vec)) ≈ vec

        vec = SA_F64[1, 1, 2]
        @test vector2spin(vec).longitudinal ≈ 2.
        @test vector2spin(vec).transverse ≈ sqrt(2)
        @test phase(vector2spin(vec)) ≈ 45.
        @test vector(vector2spin(vec)) ≈ vec
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
            snaps = evolve(Spin(), Microstructure(), Sequence(2.8), yield_every=0.5)
            time = 0.
            for snap in snaps
                @test snap.time == time
                time += 0.5
                @test vector(snap) == SA_F64[0., 0., 1.]
                @test longitudinal(snap) == 1.
                @test transverse(snap) == 0.
            end
            @test length(snaps) == 6

            snaps = evolve([Spin(), Spin()], Microstructure(), Sequence(2.8), yield_every=0.5)
            time = 0.
            for snap in snaps
                @test snap.time == time
                time += 0.5
                @test vector(snap) == SA_F64[0., 0., 2.]
                @test longitudinal(snap) == 2.
                @test transverse(snap) == 0.
            end
            @test length(snaps) == 6
        end
        @testset "Gradient echo sequence" begin
            snaps = evolve(Spin(), Microstructure(), Sequence([RFPulse(flip_angle=90)], 2.8), yield_every=0.5)
            s1 = snaps[1]
            @test vector(s1) == SA_F64[0., 0., 1.]
            time = 0.
            for snap in snaps[2:end]
                @test vector(snap) ≈ SA_F64[0., 1., 0.]
            end
            @test length(snaps) == 9
        end
        @testset "Ensure data is stored at final TR" begin
            snaps = evolve(Spin(), Microstructure(), Sequence(2.), yield_every=0.5)
            @test length(snaps) == 5
        end
        @testset "Basic diffusion has no effect in constant fields" begin
            sequence = Sequence([RFPulse(flip_angle=90)], 2.)
            no_diff = evolve(Spin(), Microstructure(R2=field(0.3)), sequence, yield_every=0.5)
            with_diff = evolve(Spin(), Microstructure(diffusivity=field(1.), R2=field(0.3)), sequence, yield_every=0.5)
            with_diff_grad = evolve(Spin(), Microstructure(diffusivity=field(1.), R2=field(0.3)), sequence, yield_every=0.5)
            spin_no_diff = no_diff[end].spins[1]
            spin_with_diff = with_diff[end].spins[1]
            @test spin_no_diff.position == SA_F64[0, 0, 0]
            @test spin_with_diff.position != SA_F64[0, 0, 0]
            @test spin_with_diff.orientation == spin_no_diff.orientation
        end
        @testset "Basic diffusion changes spin orientation in spatially varying field" begin
            sequence = Sequence([RFPulse(flip_angle=90)], 2.)
            no_diff = evolve(Spin(), Microstructure(R2=field(SA_F64[1., 0, 0], 0.3)), sequence, yield_every=0.5)
            with_diff = evolve(Spin(), Microstructure(diffusivity=field(1.), R2=field(SA_F64[1., 0., 0.], 0.3)), sequence, yield_every=0.5)
            with_diff_no_grad = evolve(Spin(), Microstructure(diffusivity=field(1.), R2=field(0.3)), sequence, yield_every=0.5)
            spin_no_diff = no_diff[end].spins[1]
            spin_with_diff = with_diff[end].spins[1]
            spin_with_diff_no_grad = with_diff_no_grad[end].spins[1]
            @test spin_no_diff.position == SA_F64[0, 0, 0]
            @test spin_with_diff.position != SA_F64[0, 0, 0]
            @test spin_with_diff_no_grad.position != SA_F64[0, 0, 0]
            @test spin_with_diff.orientation != spin_no_diff.orientation
            @test spin_with_diff_no_grad.orientation == spin_no_diff.orientation
        end
    end
end
