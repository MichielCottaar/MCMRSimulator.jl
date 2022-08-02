using Test
import MRSimulator: MRSimulator, Spin, Microstructure, evolve_to_time, time, field,
    gyromagnetic_ratio, RFPulse, apply_pulse, phase, longitudinal, transverse, time, position, 
    norm_angle, evolve, Sequence, relax, vector2spin, vector, Wall, correct_collisions, Movement,
    Cylinder, Sphere, cylinder_plane, ray_grid_intersections
using StaticArrays
using LinearAlgebra


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
    @testset "Collision tests" begin
        function compare(ms1 :: AbstractVector{Movement}, ms2 :: AbstractVector{Movement})
            @test length(ms1) == length(ms2)
            for (m1, m2) in zip(ms1, ms2)
                @test m1.origin ≈ m2.origin atol=1e-12 rtol=1e-6
                @test m1.destination ≈ m2.destination atol=1e-12 rtol=1e-6
                @test m1.timestep ≈ m2.timestep atol=1e-12 rtol=1e-6
            end
        end
        @testset "Wall reflections" begin
            @testset "Hitting vertical wall directly" begin
                res = correct_collisions(
                    Movement(SA_F64[0, 0, 0], SA_F64[3, 0, 0], 3.),
                    [Wall(:x, 1.)],
                )
                @test length(res) == 2
                @test res[1].origin ≈ SA_F64[0, 0, 0]
                @test res[1].destination ≈ SA_F64[1, 0, 0]
                @test res[1].timestep ≈ 1.
                @test res[2].origin ≈ SA_F64[1, 0, 0]
                @test res[2].destination ≈ SA_F64[-1, 0, 0]
                @test res[2].timestep ≈ 2.
            end
            @testset "Hitting vertical wall under angle" begin
                res = correct_collisions(
                    Movement(SA_F64[0, 0, 0], SA_F64[3, 6, 0], 3.),
                    [Wall(:x, 1.)],
                )
                @test length(res) == 2
                @test res[1].origin ≈ SA_F64[0, 0, 0]
                @test res[1].destination ≈ SA_F64[1, 2, 0]
                @test res[1].timestep ≈ 1.
                @test res[2].origin ≈ SA_F64[1, 2, 0]
                @test res[2].destination ≈ SA_F64[-1, 6, 0]
                @test res[2].timestep ≈ 2.
            end
            @testset "Missing vertical wall" begin
                res = correct_collisions(
                    Movement(SA_F64[0, 0, 0], SA_F64[0, 2, 0], 4.),
                    [Wall(:x, 1.)],
                )
                @test length(res) == 1
                @test res[1].origin == SA_F64[0, 0, 0]
                @test res[1].destination == SA_F64[0, 2, 0]
                @test res[1].timestep == 4.
            end
            @testset "Hitting two vertical walls" begin
                res = correct_collisions(
                    Movement(SA_F64[0, 0, 0], SA_F64[6, 0, 12], 30),
                    [Wall(:x, 1.), Wall(:x, -1.)],
                )
                @test length(res) == 4
                @test res[1].origin ≈ SA_F64[0, 0, 0]
                @test res[1].destination ≈ SA_F64[1, 0, 2]
                @test res[1].timestep ≈ 5.
                @test res[2].origin ≈ SA_F64[1, 0, 2]
                @test res[2].destination ≈ SA_F64[-1, 0, 6]
                @test res[2].timestep ≈ 10.
                @test res[3].origin ≈ SA_F64[-1, 0, 6]
                @test res[3].destination ≈ SA_F64[1, 0, 10]
                @test res[3].timestep ≈ 10.
                @test res[4].origin ≈ SA_F64[1, 0, 10]
                @test res[4].destination ≈ SA_F64[0, 0, 12]
                @test res[4].timestep ≈ 5.
            end
            @testset "Hitting vertical and horizontal walls" begin
                res = correct_collisions(
                    Movement(SA_F64[-1, 0, 0], SA_F64[2, 3, 3], 3),
                    [Wall(:x, 1.), Wall(:y, 1.)],
                )
                @test length(res) == 3
                @test res[1].origin ≈ SA_F64[-1, 0, 0]
                @test res[1].destination ≈ SA_F64[0, 1, 1]
                @test res[1].timestep ≈ 1.
                @test res[2].origin ≈ SA_F64[0, 1, 1]
                @test res[2].destination ≈ SA_F64[1, 0, 2]
                @test res[2].timestep ≈ 1.
                @test res[3].origin ≈ SA_F64[1, 0, 2]
                @test res[3].destination ≈ SA_F64[0, -1, 3]
                @test res[3].timestep ≈ 1.
            end
            @testset "Hitting a corner" begin
                res = correct_collisions(
                    Movement(SA_F64[0, 0, 0], SA_F64[3, 3, 3], 3),
                    [Wall(:x, 1.), Wall(:y, 1.)],
                )
                @test length(res) == 3
                @test res[1].origin ≈ SA_F64[0, 0, 0]
                @test res[1].destination ≈ SA_F64[1, 1, 1]
                @test res[1].timestep ≈ 1.
                @test res[2].origin ≈ SA_F64[1, 1, 1]
                @test res[2].destination ≈ SA_F64[1, 1, 1]
                @test res[2].timestep ≈ 0. atol=1e-10
                @test res[3].origin ≈ SA_F64[1, 1, 1]
                @test res[3].destination ≈ SA_F64[-1, -1, 3]
                @test res[3].timestep ≈ 2.
            end
            @testset "Hitting diagonal wall" begin
                res = correct_collisions(
                    Movement(SA_F64[1, 0, 0], SA_F64[1, 3, 0], 6.),
                    [Wall(SA_F64[1, 1, 0], 1.)],
                )
                @test length(res) == 2
                @test res[1].origin ≈ SA_F64[1, 0, 0]
                @test res[1].destination ≈ SA_F64[1, 1, 0]
                @test res[1].timestep ≈ 2.
                @test res[2].origin ≈ SA_F64[1, 1, 0]
                @test res[2].destination ≈ SA_F64[-1, 1, 0]
                @test res[2].timestep ≈ 4.
            end
        end
        @testset "Sphere reflections" begin
            @testset "Remain within sphere" begin
                res = correct_collisions(
                    Movement(SA_F64[0, 0, 1.5], SA_F64[6, 8, 6], 6.),
                    [Sphere(1., SA_F64[0, 0, 2])]
                )
                final = res[end].destination
                radius = norm(final .- SA_F64[0, 0, 2])
                @assert radius <= 1.
            end
        end
        @testset "Cylinder reflections" begin
            @testset "Within cylinder along radial line" begin
                res = correct_collisions(
                    Movement(SA_F64[0, 0, 0], SA_F64[6, 0, 6], 6.),
                    [Cylinder(1., :z, SA_F64[0, 0, 2])]
                )
                @test length(res) == 4
                @test res[1].origin ≈ SA_F64[0, 0, 0]
                @test res[1].destination ≈ SA_F64[1, 0, 1]
                @test res[1].timestep ≈ 1.
                @test res[2].origin ≈ SA_F64[1, 0, 1]
                @test res[2].destination ≈ SA_F64[-1, 0, 3]
                @test res[2].timestep ≈ 2.
                @test res[3].origin ≈ SA_F64[-1, 0, 3]
                @test res[3].destination ≈ SA_F64[1, 0, 5]
                @test res[3].timestep ≈ 2.
                @test res[4].origin ≈ SA_F64[1, 0, 5]
                @test res[4].destination ≈ SA_F64[0, 0, 6]
                @test res[4].timestep ≈ 1.
            end
            @testset "90 degree bounces within vertical cylinder" begin
                res = correct_collisions(
                    Movement(SA_F64[0, 1, 0], SA_F64[10, 1, 0], 10),
                    [Cylinder(sqrt(2), :z, SA_F64[0, 0, 2])]
                )
                compare(res, [
                    Movement(SA_F64[0, 1, 0], SA_F64[1, 1, 0], 1),
                    Movement(SA_F64[1, 1, 0], SA_F64[1, -1, 0], 2),
                    Movement(SA_F64[1, -1, 0], SA_F64[-1, -1, 0], 2),
                    Movement(SA_F64[-1, -1, 0], SA_F64[-1, 1, 0], 2),
                    Movement(SA_F64[-1, 1, 0], SA_F64[1, 1, 0], 2),
                    Movement(SA_F64[1, 1, 0], SA_F64[1, 0, 0], 1),
                ])
            end
            @testset "90 degree bounces from outside vertical cylinder" begin
                res = correct_collisions(
                    Movement(SA_F64[-2, 1, 0], SA_F64[2, 1, 0], 4),
                    [Cylinder(sqrt(2), :z, SA_F64[0, 0, 2])]
                )
                compare(res, [
                    Movement(SA_F64[-2, 1, 0], SA_F64[-1, 1, 0], 1),
                    Movement(SA_F64[-1, 1, 0], SA_F64[-1, 4, 0], 3),
                ])
            end
            @testset "Remain within angled cylinder" begin
                orient = SA_F64[1, 2, sqrt(3)]
                cylinder = Cylinder(2.3, orient, SA_F64[0, 0, 0])
                res = correct_collisions(
                    Movement(SA_F64[0, 0.5, 0.3], SA_F64[-30, 50, 10], 40),
                    [cylinder]
                )
                final = res[end].destination
                radius = norm(final .- (orient ⋅ final) * orient / norm(orient) ^ 2)
                @test radius <= 2.3
            end
        end
        @testset "Ray-grid intersections" begin
            function tcompare(t1, t2)
                @test length(t1) == length(t2)
                for (e1, e2) in zip(t1, t2)
                    @test e1 ≈ e2
                end
            end
            res = collect(MRSimulator.ray_grid_intersections(SA_F64[0.5, 0.5, 0.5], SA_F64[0.5, 0.5, 3.5]))
            tcompare(res[1], ([0, 0, 0], 0., [0.5, 0.5, 0.5], 1/6, [0.5, 0.5, 1.]))
            tcompare(res[2], ([0, 0, 1], 1/6, [0.5, 0.5, 0.], 1/2, [0.5, 0.5, 1.]))
            tcompare(res[3], ([0, 0, 2], 1/2, [0.5, 0.5, 0.], 5/6, [0.5, 0.5, 1.]))
            tcompare(res[4], ([0, 0, 3], 5/6, [0.5, 0.5, 0.], 1., [0.5, 0.5, 0.5]))
        end
        @testset "Reflections on planes of cylinders" begin
            @testset "Bounce between four cylinders" begin
                cylinders = cylinder_plane(
                    sqrt(2), repeatx=3, repeaty=4
                )
                res = correct_collisions(
                    Movement(SA_F64[1, 2, 0], SA_F64[1, 11, 9], 9),
                    [cylinders]
                )
                compare(res, [
                    Movement(SA_F64[1, 2, 0], SA_F64[1, 3, 1], 1),
                    Movement(SA_F64[1, 3, 1], SA_F64[2, 3, 2], 1),
                    Movement(SA_F64[2, 3, 2], SA_F64[2, 1, 4], 2),
                    Movement(SA_F64[2, 1, 4], SA_F64[1, 1, 5], 1),
                    Movement(SA_F64[1, 1, 5], SA_F64[1, 3, 7], 2),
                    Movement(SA_F64[1, 3, 7], SA_F64[2, 3, 8], 1),
                    Movement(SA_F64[2, 3, 8], SA_F64[2, 2, 9], 1),
                ])
            end
            @testset "Travel through many repeats between bounces" begin
                cylinders = cylinder_plane(
                    1, repeatx=2, repeaty=4
                )
                res = correct_collisions(
                    Movement(SA_F64[1, 2, 0], SA_F64[13, 6, 4], 12),
                    [cylinders]
                )
                compare(res, [
                    Movement(SA_F64[1, 2, 0], SA_F64[4, 3, 1], 3),
                    Movement(SA_F64[4, 3, 1], SA_F64[10, 1, 3], 6),
                    Movement(SA_F64[10, 1, 3], SA_F64[13, 2, 4], 3),
                ])
            end
        end
    end
end
