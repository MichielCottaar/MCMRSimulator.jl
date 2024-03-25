@testset "test_evolve.jl" begin
    @testset "Timings of readouts" begin
        @testset "Explicitly setting timestep" begin
            s1 = SpinEcho(TE=80.)

            eval_at = [
                0., nextfloat(0.),
                0.5:0.5:39.5...,
                prevfloat(40.), 40., nextfloat(40.),
                40.5:0.5:79.5...,
                prevfloat(80.), 80.,
            ]
            @test all(mr.propose_times(mr.Simulation(s1, max_timestep=0.5), 0., 80.) .≈ eval_at)
            @test all(mr.propose_times(mr.Simulation(s1, max_timestep=0.501), 0., 80.) .≈ eval_at)
        end
        @testset "Setting gradient_precision" begin
            s1 = SpinEcho(TE=80.)
            s2 = DWI(TE=80., bval=1.)

            kwargs = Dict(
                :diffusivity => 0.,
                :gradient_precision => 1.,
                :max_timestep => Inf,
            )
            control_points = [0., nextfloat(0.), prevfloat(40.), 40., nextfloat(40.), prevfloat(80.), 80.]

            @test all(mr.propose_times(mr.Simulation(s1; kwargs...), 0., 80.) .== control_points)
            @test all(mr.propose_times(mr.Simulation(s1; kwargs...), 90., 120.) .== [90., 120.])
            @test all(mr.propose_times(mr.Simulation([s1, s2]; kwargs...), 0., 80.) .== control_points)
            @test all(mr.propose_times(mr.Simulation([s1, s2]; kwargs...), 90., 120.) .== [90., 120.])

            kwargs[:diffusivity] = 1.
            @test all(mr.propose_times(mr.Simulation(s1; kwargs...), 0., 80.) .== control_points)
            @test all(mr.propose_times(mr.Simulation(s1; kwargs...), 90., 120.) .== [90., 120.])
            @test length(mr.propose_times(mr.Simulation([s1, s2]; kwargs...), 0., 80.)) > length(control_points)
            @test length(intersect(mr.propose_times(mr.Simulation([s1, s2]; kwargs...), 0., 80.), control_points)) == length(control_points)
            @test all(mr.propose_times(mr.Simulation([s1, s2]; kwargs...), 90., 120.) .== [90., 120.])
        end
    end
    @testset "Empty environment and sequence" begin
        empty_sequence = build_sequence() do
            Sequence([2.8])
        end
        simulation = mr.Simulation(empty_sequence)
        snaps = mr.readout(zeros(3), simulation, 0:0.5:2.8, return_snapshot=true)
        time = 0.
        for snap in snaps
            @test snap.time == time
            time += 0.5
            @test mr.orientation(snap) == SA[0., 0., 1.]
            @test mr.longitudinal(snap) == 1.
            @test mr.transverse(snap) == 0.
        end
        @test length(snaps) == 6

        simulation = mr.Simulation(empty_sequence)
        snaps= mr.readout([mr.Spin(), mr.Spin()], simulation, 0:0.5:2.8, return_snapshot=true)
        time = 0.
        for snap in snaps
            @test snap.time == time
            time += 0.5
            @test mr.orientation(snap) == SA[0., 0., 2.]
            @test mr.longitudinal(snap) == 2.
            @test mr.transverse(snap) == 0.
        end
        @test length(snaps) == 6
    end
    @testset "Gradient echo sequence" begin
        simulation = mr.Simulation(GradientEcho(TE=2.8))
        snaps = mr.readout(zeros(3), simulation, 0:0.5:2.8)
        @test mr.orientation(snaps[1]) ≈ SA[0., 0., 1.]
        for snap in snaps[2:end]
            @test mr.orientation(snap) ≈ SA[0., 1., 0.]
        end
        @test length(snaps) == 6
    end
    @testset "Ensure data is stored at requested time" begin
        empty_sequence = build_sequence() do
            Sequence([2.8])
        end
        simulation = mr.Simulation(empty_sequence)

        snaps = mr.evolve(mr.Spin(), simulation, 2.3)
        @test mr.get_time(snaps) == 2.3

        snaps = mr.evolve(snaps, simulation)
        @test mr.get_time(snaps) == 2.8

        snaps = mr.evolve(snaps, simulation)
        @test mr.get_time(snaps) == 5.6

        snaps = mr.evolve(mr.Spin(), simulation)
        @test mr.get_time(snaps) == 2.8
    end
    @testset "Basic diffusion has no effect in constant fields" begin
        sequence = GradientEcho(TE=2.)
        no_diff = mr.Simulation([sequence], diffusivity=0., R2=0.3)
        with_diff = mr.Simulation([sequence], diffusivity=1., R2=0.3)
        spin_no_diff = mr.evolve(mr.Spin(), no_diff).spins[1]
        spin_with_diff = mr.evolve(mr.Spin(), with_diff).spins[1]
        @test spin_no_diff.position == SA[0, 0, 0]
        @test spin_with_diff.position != SA[0, 0, 0]
        @test mr.orientation.(spin_with_diff.orientations) == mr.orientation.(spin_no_diff.orientations)
        @test mr.transverse(spin_no_diff) ≈ exp(-0.6)
        @test abs(mr.longitudinal(spin_no_diff)) < Float64(1e-6)
    end
    @testset "Basic diffusion run within sphere" begin
        sequence = GradientEcho(TE=20.)
        sphere = mr.Spheres(radius=1.)
        Random.seed!(12)
        diff = mr.Simulation(sequence, diffusivity=2., geometry=sphere)
        snaps = mr.readout([mr.Spin(), mr.Spin()], diff, 0:0.5:sequence.TR, return_snapshot=true)
        @test size(snaps) == (5, )
        for snap in snaps
            @test length(snap.spins) == 2
            for spin in snap.spins
                @test norm(spin.position) < 1.
                @test length(spin.orientations) == 1
            end
        end
    end
    @testset "Run simulation with multiple sequences at once" begin
        sequences = [
            build_sequence() do 
                Sequence([InstantPulse(flip_angle=0, phase=0.), 2., SingleReadout(), 1.]) 
            end,
            build_sequence() do 
                Sequence([InstantPulse(flip_angle=90, phase=0.), 2., SingleReadout(), 1.]) 
            end,
            build_sequence() do 
                Sequence([InstantPulse(flip_angle=90, phase=0.), 1., SingleReadout(), 1.]) 
            end,
        ]
        all_snaps = mr.Simulation(sequences, diffusivity=1., R2=1.)

        readouts = mr.readout(mr.Spin(), all_snaps)
        @assert size(readouts) == (3,)

        # check relaxation
        @test mr.transverse(readouts[1]) ≈ 0. atol=1e-12
        @test mr.transverse(readouts[2]) ≈ exp(-2.)
        @test mr.transverse(readouts[3]) ≈ exp(-1.)
        @test mr.longitudinal(readouts[1]) ≈ 1.
        @test mr.longitudinal(readouts[2]) ≈ 0. atol=1e-12
        @test mr.longitudinal(readouts[3]) ≈ 0. atol=1e-12
    end
end
