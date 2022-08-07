using Test
import MRSimulator: MRSimulator, Spin, Microstructure, evolve_to_time, time, field,
    gyromagnetic_ratio, RFPulse, apply, phase, longitudinal, transverse, time, position, 
    norm_angle, evolve_TR, Sequence, relax, vector2spin, vector, Wall, correct_collisions, Movement,
    Cylinder, Sphere, cylinder_plane, ray_grid_intersections, StepTrajectory, Obstruction, InstantGradient
import MRSimulator as mr
using StaticArrays
using LinearAlgebra
import Random


@testset "MRSimulator tests" begin
    include("test_field.jl")
    include("test_collisions.jl")
    include("test_evolve.jl")
    include("test_known_sequences.jl")
    @testset "Simple relaxation" begin
        orient = Spin(transverse=1., longitudinal=0.).orientation
        pos = zero(SVector{3, Float64})
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
    @testset "Apply Sequence components" begin
        @testset "0 degree pulses should do nothing" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = RFPulse(0., 0., pulse_phase)
                for spin_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                    spin = Spin(phase=spin_phase, transverse=1.)
                    spin = apply(pulse, spin.orientation)
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
                spin = apply(pulse, spin.orientation)
                @test longitudinal(spin) ≈ -1.
            end
        end
        @testset "90 degree pulses should eliminate longitudinal" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = RFPulse(0., 90., pulse_phase)
                spin = Spin()
                @test longitudinal(spin) == 1.
                spin = apply(pulse, spin.orientation)
                @test longitudinal(spin) ≈ 0. atol=1e-12
            end
        end
        @testset "Spins with same phase as pulse are unaffected by pulse" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                for flip_angle in (10, 90, 120, 180)
                    pulse = RFPulse(0., flip_angle, pulse_phase)

                    spin = Spin(longitudinal=0., transverse=1., phase=pulse_phase)
                    spin = apply(pulse, spin.orientation)
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
                spin = apply(pulse, spin.orientation)
                @test longitudinal(spin) ≈ 0. atol=1e-12
                @test transverse(spin) ≈ 1.
                @test phase(spin) ≈ -spin_phase
            end
            for pulse_phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 180, pulse_phase)

                spin = Spin(longitudinal=0., transverse=1., phase=0.)
                spin = apply(pulse, spin.orientation)
                @test longitudinal(spin) ≈ 0. atol=1e-12
                @test transverse(spin) ≈ 1.
                @test phase(spin) ≈ 2 * pulse_phase
            end
        end
        @testset "90 pulses flips longitudinal spin into transverse plane" begin
            for pulse_phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 90, pulse_phase)
                spin = apply(pulse, Spin().orientation)
                @test longitudinal(spin) ≈ 0. atol=1e-12
                @test transverse(spin) ≈ 1.
                @test phase(spin) ≈ pulse_phase + 90
            end
        end
        @testset "90 pulses flips transverse spin into longitudinal plane" begin
            for pulse_phase in (0., 22., 30., 80.)
                pulse = RFPulse(0., 90, pulse_phase)
                spin_phase = (pulse_phase + 90)
                spin = Spin(longitudinal=0., transverse=1., phase=spin_phase)
                spin = apply(pulse, spin.orientation)
                @test longitudinal(spin) ≈ -1.
                @test transverse(spin) ≈ 0. atol=1e-12
            end
        end
        @testset "Gradient should do nothing at origin" begin
            spin = Spin(position=SA_F64[0, 2, 2], transverse=1., phase=90.)
            @test phase(spin) == 90.
            spin2 = apply(InstantGradient(qvec=SA_F64[4, 0, 0]), spin)
            @test phase(spin2) == 90.
        end
        @testset "Test instant gradient effect away from origin" begin
            spin = Spin(position=SA_F64[2, 2, 2], transverse=1., phase=90.)
            @test phase(spin) == 90.
            spin2 = apply(InstantGradient(qvec=SA_F64[0.25, 0, 0]), spin)
            @test phase(spin2) == 90. + rad2deg(0.5)
        end
    end
    @testset "Random generator number control" begin
        @testset "FixedXoshiro interface" begin
            rng = mr.FixedXoshiro()
            copy!(Random.TaskLocalRNG(), rng)
            a = randn(2)
            copy!(Random.TaskLocalRNG(), rng)
            @test randn(2) == a
        end
    end
end
