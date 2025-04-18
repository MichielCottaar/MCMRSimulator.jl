
@testset "test_various.jl" begin
@test length(detect_ambiguities(mr)) == 0

@testset "Draw random positions and radii" begin
    pos, rad = mr.random_positions_radii([20., 20.], 0.7, 2; variance=0.01)
    n_too_close = 0
    for i in 1:length(pos)
        for j in i+1:length(pos)
            dist = norm(@. mod(pos[i] - pos[j] + 15, 30) - 15)
            if dist < (rad[i] + rad[j])
                n_too_close += 1
            end
        end
    end
    @test n_too_close == 0
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
            pulse = InstantPulse(0., 0., nothing)
            for spin_phase in (-90, -45, 0., 30., 90., 170, 22.123789)
                spin = mr.Spin(phase=spin_phase, transverse=1.)
                mr.Evolve.apply_instants!([spin], 1, pulse, nothing)
                @test mr.phase(spin) ≈ spin_phase
                @test mr.longitudinal(spin) ≈ 1.
                @test mr.transverse(spin) ≈ 1.
            end
        end
    end
    @testset "180 degree pulses should flip longitudinal" begin
        for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
            pulse = InstantPulse(180., pulse_phase, nothing)
            spin = mr.Spin()
            @test mr.longitudinal(spin) == 1.
            mr.Evolve.apply_instants!([spin], 1, pulse, nothing)
            @test mr.longitudinal(spin) ≈ -1.
        end
    end
    @testset "90 degree pulses should eliminate longitudinal" begin
        for pulse_phase in (-90, -45, 0., 30., 90., 180, 22.123789)
            pulse = InstantPulse(90., pulse_phase, nothing)
            spin = mr.Spin()
            @test mr.longitudinal(spin) == 1.
            mr.Evolve.apply_instants!([spin], 1, pulse, nothing)
            @test mr.longitudinal(spin) ≈ 0. atol=1e-7
        end
    end
    @testset "Spins with same phase as pulse are unaffected by pulse" begin
        for pulse_phase in (-90, -45, 0., 30., 90., 170, 22.123789)
            for flip_angle in (10, 90, 120, 180)
                pulse = InstantPulse(flip_angle, pulse_phase, nothing)

                spin = mr.Spin(longitudinal=0., transverse=1., phase=pulse_phase)
                mr.Evolve.apply_instants!([spin], 1, pulse, nothing)
                @test mr.longitudinal(spin) ≈ zero(Float64) atol=1e-7
                @test mr.transverse(spin) ≈ one(Float64)
                @test mr.phase(spin) ≈ Float64(pulse_phase)
            end
        end
    end
    @testset "180 pulses flips phase around axis" begin
        for spin_phase in (0., 22., 30., 80.)
            pulse = InstantPulse(180, 0., nothing)

            spin = mr.Spin(longitudinal=0., transverse=1., phase=spin_phase)
            mr.Evolve.apply_instants!([spin], 1, pulse, nothing)
            @test mr.longitudinal(spin) ≈ 0. atol=1e-7
            @test mr.transverse(spin) ≈ 1.
            @test mr.phase(spin) ≈ -spin_phase
        end
        for pulse_phase in (0., 22., 30., 80.)
            pulse = InstantPulse(180, pulse_phase, nothing)

            spin = mr.Spin(longitudinal=0., transverse=1., phase=0.)
            mr.Evolve.apply_instants!([spin], 1, pulse, nothing)
            @test mr.longitudinal(spin) ≈ 0. atol=1e-7
            @test mr.transverse(spin) ≈ 1.
            @test mr.phase(spin) ≈ 2 * pulse_phase
        end
    end
    @testset "90 pulses flips longitudinal spin into transverse plane" begin
        for pulse_phase in (0., 22., 30., 80.)
            pulse = InstantPulse(90, pulse_phase, nothing)
            spin = mr.Spin()
            mr.Evolve.apply_instants!([spin], 1, pulse, nothing)
            @test mr.longitudinal(spin) ≈ 0. atol=1e-7
            @test mr.transverse(spin) ≈ 1.
            @test mr.phase(spin) ≈ pulse_phase - 90
        end
    end
    @testset "90 pulses flips transverse spin into longitudinal plane" begin
        for pulse_phase in (0., 22., 30., 80.)
            pulse = InstantPulse(90, pulse_phase, nothing)
            spin_phase = (pulse_phase - 90)
            spin = mr.Spin(longitudinal=0., transverse=1., phase=spin_phase)
            mr.Evolve.apply_instants!([spin], 1, pulse, nothing)
            @test mr.longitudinal(spin) ≈ -1.
            @test mr.transverse(spin) ≈ 0. atol=1e-7
        end
    end
    @testset "Gradient should do nothing at origin" begin
        spin = mr.Spin(position=[0, 2, 2], transverse=1., phase=90.)
        @test mr.phase(spin) ≈ Float64(90.)
        mr.Evolve.apply_instants!([spin], 1, InstantGradient(qvec=[4, 0, 0]), nothing)
        @test mr.phase(spin) ≈ Float64(90.)
    end
    @testset "Test instant gradient effect away from origin" begin
        spin = mr.Spin(position=[2, 2, 2], transverse=1., phase=90.)
        @test mr.phase(spin) ≈ Float64(90.)
        mr.Evolve.apply_instants!([spin], 1, InstantGradient(qvec=[0.01, 0, 0]), nothing)
        @test mr.phase(spin) ≈ Float64(90. + 0.02 * 360 / 2π)
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
        t1 = mr.readout(spin, env, 1:5, return_snapshot=true)
        t2 = mr.readout(spin, env, 1:5, return_snapshot=true)
        @test all(@. mr.position(t1) == mr.position(t2))
        get_rng(snapshot) = snapshot.spins[1].rng
        @test all(@. get_rng(t1) == get_rng(t2))
    end
end
@testset "Test readout formats" begin
    sequence = GradientEcho(TE=1000)
    seq = build_sequence() do 
        Sequence([10., SingleReadout(), 20., SingleReadout(), 70]) 
    end
    sim_empty = mr.Simulation([])
    sim_flat = mr.Simulation(seq)
    sim_single = mr.Simulation([seq])
    sim_double = mr.Simulation([seq, seq])

    @testset "Test readout function output" begin
        times = 1:0.1:10
        Nt = length(times)
        for (sim, shape) in [
            (sim_flat, ()),
            (sim_single, (1, )),
            (sim_double, (2, )),
        ]
            @test size(mr.readout(10, sim, times, return_snapshot=true)) == (shape..., Nt)
            @test size(mr.readout(10, sim, times, return_snapshot=true, nTR=1)) == (shape..., Nt, 1)
            @test size(mr.readout(10, sim, times, return_snapshot=true, nTR=2)) == (shape..., Nt, 2)
            @test size(mr.readout(10, sim, return_snapshot=true)) == (shape..., 2)
            @test size(mr.readout(10, sim, nTR=1)) == (shape..., 2, 1)
            @test size(mr.readout(10, sim, nTR=2)) == (shape..., 2, 2)
            @test all(mr.longitudinal.(mr.readout(10, sim, times)) .≈ 10.)
        end

        @test_throws ErrorException mr.readout(100, sim_empty)
        @test_throws ErrorException mr.readout(100, sim_empty, return_snapshot=true)
        @test_throws ErrorException mr.readout(100, sim_empty, times)
        @test size(mr.readout(100, sim_empty, times, return_snapshot=true)) == (Nt, )
        @test_throws ErrorException mr.readout(100, sim_empty, times, return_snapshot=true, nTR=1)
        @test_throws ErrorException mr.readout(100, sim_empty, times, return_snapshot=true, nTR=2)

        @test size(mr.readout(100, sim_flat)) == (2, )
        for (time, snapshot) in zip((10., 30.), mr.readout(100, sim_flat, return_snapshot=true))
            @test length(snapshot) == 100
            @test time == mr.get_time(snapshot)
            @test mr.longitudinal(snapshot) ≈ 100.
        end
        @test size(mr.readout(100, sim_single)) == (1, 2)
        @test size(mr.readout(100, sim_double)) == (2, 2)
        for sim in (sim_single, sim_double)
            @test length(sim.sequences) == size(mr.readout(100, sim))[1]
            for r in eachrow(mr.readout(100, sim, return_snapshot=true))
                @test length(r) == 2
                for (time, snapshot) in zip((10., 30.), r)
                    @test length(snapshot) == 100
                    @test time == mr.get_time(snapshot)
                    @test mr.longitudinal(snapshot) ≈ 100.
                end
            end
        end
    end
end

@testset "Test simulation pretty printing" begin
    sim = mr.Simulation([DWI(TE=80., bval=1), DWI(TE=80., bval=2)], geometry=mr.Spheres(radius=[1, 2.], repeats=[5, 5, 5]), R1=0.1, surface_relaxation=0.3)
    @test repr(sim, context=:compact => true) == "Simulation(2 sequences, Geometry(2 repeating Round objects, ), D=3.0um^2/ms, GlobalProperties(R1=0.1kHz, ))" 
end

@testset "Test size scale calculations" begin
    size_scale(g) = mr.Geometries.Internal.size_scale(mr.fix(g))

    @test isinf(size_scale(mr.Walls(position=0)))
    @test size_scale(mr.Walls(position=[0, 2])) == 2.
    @test size_scale(mr.Walls(repeats=5)) == 5.
    @test size_scale(mr.Walls(position=[0, 2], repeats=5)) == 2.
    @test size_scale(mr.Walls(position=[0, 4], repeats=5)) == 1.

    @test size_scale(mr.Cylinders(radius=[0.3, 0.8], repeats=[2, 3])) == 0.3
    @test size_scale(mr.Spheres(radius=[0.3, 0.8])) == 0.3
    @test size_scale(mr.Annuli(inner=[0.5, 0.7], outer=[0.6, 0.9])) ≈ 0.1
end

@testset "Test geometry JSON I/O" begin
    for geometry in (
        mr.Cylinders(radius=[0.7, 0.8]),
        mr.Annuli(inner=0.2, outer=0.4, repeats=[2, 3], rotation=:y),
        mr.Walls(position=[1, 2, 3, 4, 5]),
        mr.Spheres(radius=[0.2, 0.4, 0.8], R1_inside=1., R2_inside=[2., 2.3, 3.4]),
        mr.Mesh(triangles=[[1, 2, 3], [2, 1, 3]], vertices=[[0, 0, 0], [1, 1, 1], [2, 0, 3]])
    )
        io = IOBuffer()
        mr.write_geometry(io, geometry)
        s = String(take!(io))
        group = mr.read_geometry(s)
        @test group.n_obstructions == geometry.n_obstructions
        for (key, field_value) in group.field_values
            @test all(field_value.value .== getproperty(geometry, key).value)
        end

        io = IOBuffer()
        mr.write_geometry(io, [geometry])
        s = String(take!(io))
        groups = mr.read_geometry(s)
        @test length(groups) == 1
        group = groups[1]
        @test group.n_obstructions == geometry.n_obstructions
        for (key, field_value) in group.field_values
            @test all(field_value.value .== getproperty(geometry, key).value)
        end
    end
end

@testset "Test splitting RF pulse application" begin
    # relevant for collisions
    props = mr.Geometries.Internal.MRIProperties(0., 0., 0.)
    cp = mr.SequenceParts.ConstantPulse(1., 0., 2.)

    orient_single = mr.SpinOrientation([0., 0., 1.])
    orient_dual = mr.SpinOrientation([0., 0., 1.])
    mr.Relax.apply_pulse!(orient_single, cp, props, 0.1, 0., 1., 0., Val(1))

    mr.Relax.apply_pulse!(orient_dual, cp, props, 0.1, 0., 0.3, 0., Val(1))
    mr.Relax.apply_pulse!(orient_dual, cp, props, 0.1, 0.3, 1., 0., Val(1))
    @test orient_single.longitudinal ≈ orient_dual.longitudinal
    @test orient_single.transverse ≈ orient_dual.transverse
    @test orient_single.phase ≈ orient_dual.phase

    orient_many = mr.SpinOrientation([0., 0., 1.])
    mr.Relax.apply_pulse!(orient_many, cp, props, 0.1, 0., 1., 0., Val(10))
    @test orient_single.longitudinal ≈ orient_many.longitudinal
    @test orient_single.transverse ≈ orient_many.transverse
    @test orient_single.phase ≈ orient_many.phase
end

end