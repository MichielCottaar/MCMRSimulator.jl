using Test
import MCMRSimulator as mr
import MCMRSimulator: Float, SA
using StaticArrays
using LinearAlgebra
import Random
import SpecialFunctions: erf
using Statistics

all_tests = [
    "collisions",
    "evolve",
    "plots",
    "field",
    "known_sequences",
    "meshes",
    "offresonance",
    "transfer",
    "permeability",
]

if length(ARGS) == 0
    tests = all_tests
elseif length(ARGS) == 1 && ARGS[1] == "no-plots"
    tests = symdiff(all_tests, ["plots"])
else
    tests = ARGS
end


@testset "MCMRSimulator tests" begin
    for test in tests
        if test == "plots"
            include("visual_tests/run_visual_tests.jl")
        else
            include("test_$test.jl")
        end
    end
    @test length(detect_ambiguities(mr)) == 0
    @testset "Simple relaxation" begin
        orient = mr.Spin(transverse=1., longitudinal=0.).orientations[1]
        pos = zero(SVector{3, Float})
        @testset "R2 relaxation" begin
            env = mr.Microstructure(R2=mr.field(2.))(pos)
            @test mr.relax(orient, env, 0.3).transverse ≈ exp(-0.6)
        end
        @testset "R1 relaxation" begin
            env = mr.Microstructure(R1=mr.field(2.))(pos)
            @test mr.relax(orient, env, 0.3).longitudinal ≈ 1 - exp(-0.6)
        end
    end
    @testset "Spin conversions" begin
        vec = SA[1, 0, 0]
        @test mr.SpinOrientation(vec).longitudinal ≈ 0.
        @test mr.SpinOrientation(vec).transverse ≈ 1.
        @test mr.phase(mr.SpinOrientation(vec)) ≈ 0.
        @test mr.orientation(mr.SpinOrientation(vec)) ≈ vec

        vec = SA[0, 1, 1]
        @test mr.SpinOrientation(vec).longitudinal ≈ 1.
        @test mr.SpinOrientation(vec).transverse ≈ 1.
        @test mr.phase(mr.SpinOrientation(vec)) ≈ 90.
        @test mr.orientation(mr.SpinOrientation(vec)) ≈ vec

        vec = SA[1, 1, 2]
        @test mr.SpinOrientation(vec).longitudinal ≈ 2.
        @test mr.SpinOrientation(vec).transverse ≈ sqrt(2)
        @test mr.phase(mr.SpinOrientation(vec)) ≈ 45.
        @test mr.orientation(mr.SpinOrientation(vec)) ≈ vec
    end
    @testset "Apply Sequence components" begin
        @testset "0 degree flip angle pulses should do nothing" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = mr.RFPulse(0., 0., pulse_phase)
                for spin_phase in (-90, -45, 0., 30., 90., 170, 22.123789)
                    spin = mr.Spin(phase=spin_phase, transverse=1.)
                    spin = mr.apply(pulse, spin)
                    @test mr.phase(spin) ≈ spin_phase
                    @test mr.longitudinal(spin) ≈ 1.
                    @test mr.transverse(spin) ≈ 1.
                end
            end
        end
        @testset "180 degree pulses should flip longitudinal" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = mr.RFPulse(0., 180., pulse_phase)
                spin = mr.Spin()
                @test mr.longitudinal(spin) == 1.
                spin = mr.apply(pulse, spin)
                @test mr.longitudinal(spin) ≈ -1.
            end
        end
        @testset "90 degree pulses should eliminate longitudinal" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
                pulse = mr.RFPulse(0., 90., pulse_phase)
                spin = mr.Spin()
                @test mr.longitudinal(spin) == 1.
                spin = mr.apply(pulse, spin)
                @test mr.longitudinal(spin) ≈ 0. atol=1e-7
            end
        end
        @testset "Spins with same phase as pulse are unaffected by pulse" begin
            for pulse_phase in (-90, -45, 0., 30., 90., 170, 22.123789)
                for flip_angle in (10, 90, 120, 180)
                    pulse = mr.RFPulse(0., flip_angle, pulse_phase)

                    spin = mr.Spin(longitudinal=0., transverse=1., phase=pulse_phase)
                    spin = mr.apply(pulse, spin)
                    @test mr.longitudinal(spin) ≈ zero(Float) atol=1e-7
                    @test mr.transverse(spin) ≈ one(Float)
                    @test mr.phase(spin) ≈ Float(pulse_phase)
                end
            end
        end
        @testset "180 pulses flips phase around axis" begin
            for spin_phase in (0., 22., 30., 80.)
                pulse = mr.RFPulse(0., 180, 0.)

                spin = mr.Spin(longitudinal=0., transverse=1., phase=spin_phase)
                spin = mr.apply(pulse, spin)
                @test mr.longitudinal(spin) ≈ 0. atol=1e-7
                @test mr.transverse(spin) ≈ 1.
                @test mr.phase(spin) ≈ -spin_phase
            end
            for pulse_phase in (0., 22., 30., 80.)
                pulse = mr.RFPulse(0., 180, pulse_phase)

                spin = mr.Spin(longitudinal=0., transverse=1., phase=0.)
                spin = mr.apply(pulse, spin)
                @test mr.longitudinal(spin) ≈ 0. atol=1e-7
                @test mr.transverse(spin) ≈ 1.
                @test mr.phase(spin) ≈ 2 * pulse_phase
            end
        end
        @testset "90 pulses flips longitudinal spin into transverse plane" begin
            for pulse_phase in (0., 22., 30., 80.)
                pulse = mr.RFPulse(0., 90, pulse_phase)
                spin = mr.apply(pulse, mr.Spin())
                @test mr.longitudinal(spin) ≈ 0. atol=1e-7
                @test mr.transverse(spin) ≈ 1.
                @test mr.phase(spin) ≈ pulse_phase + 90
            end
        end
        @testset "90 pulses flips transverse spin into longitudinal plane" begin
            for pulse_phase in (0., 22., 30., 80.)
                pulse = mr.RFPulse(0., 90, pulse_phase)
                spin_phase = (pulse_phase + 90)
                spin = mr.Spin(longitudinal=0., transverse=1., phase=spin_phase)
                spin = mr.apply(pulse, spin)
                @test mr.longitudinal(spin) ≈ -1.
                @test mr.transverse(spin) ≈ 0. atol=1e-7
            end
        end
        @testset "Gradient should do nothing at origin" begin
            spin = mr.Spin(position=SA[0, 2, 2], transverse=1., phase=90.)
            @test mr.phase(spin) ≈ Float(90.)
            spin2 = mr.apply(mr.InstantGradient(qvec=SA[4, 0, 0]), spin)
            @test mr.phase(spin2) ≈ Float(90.)
        end
        @testset "Test instant gradient effect away from origin" begin
            spin = mr.Spin(position=SA[2, 2, 2], transverse=1., phase=90.)
            @test mr.phase(spin) ≈ Float(90.)
            spin2 = mr.apply(mr.InstantGradient(qvec=SA[0.25, 0, 0]), spin)
            @test mr.phase(spin2) ≈ Float(90. + mr.rad2deg(0.5))
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
        @testset "FixedXoshiro predictability" begin
            Random.seed!(1234)
            rng = mr.FixedXoshiro()
            Random.seed!(1234)
            rng2 = mr.FixedXoshiro()
            @test rng == rng2
        end
        @testset "Reproducible evolution" begin
            spin = mr.Spin()
            env = mr.Simulation([], diffusivity=3.)
            t1 = [s[1] for s in mr.trajectory(spin, env, 1:5)]
            t2 = [s[1] for s in mr.trajectory(spin, env, 1:5)]
            @test all(t1 .== t2)
        end
    end
    @testset "Bounding boxes" begin
        # single obstruction
        @test mr.BoundingBox(mr.Cylinder(1)) == mr.BoundingBox([-1, -1], [1, 1])

        # shifted cylinder
        @test mr.BoundingBox(mr.cylinders(1., positions=[2., 2.])) == mr.BoundingBox([1., 1.], [3, 3])

        # repeated obstructions
        @test mr.BoundingBox(mr.cylinders(1., repeats=[1., 3.])) == mr.BoundingBox([-1, -1], [1, 1])

        # shifted spheres
        @test mr.BoundingBox(mr.spheres(1., positions=[[1, 0, 0], [0, 1, 0]])) == mr.BoundingBox([-1, -1, -1.], [2., 2., 1.])
    end
end
