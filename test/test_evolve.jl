@testset "test_evolve.jl" begin
    @testset "Empty environment and sequence" begin
        empty_sequence = build_sequence(2.8)
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
        simulation = mr.Simulation(mr.read_pulseq(joinpath(@__DIR__, "pulseq", "gradient_echo_TE_2.8.seq")))
        snaps = mr.readout(zeros(3), simulation, 0:0.5:2.8)
        @test mr.orientation(snaps[1]) ≈ SA[0., 0., 1.]
        for snap in snaps[2:end]
            @test mr.orientation(snap) ≈ SA[0., -1., 0.]
        end
        @test length(snaps) == 6
    end
    @testset "Ensure data is stored at requested time" begin
        empty_sequence = build_sequence(2.8)
        simulation = mr.Simulation(empty_sequence)

        snaps = mr.evolve(mr.Spin(), simulation, 2.3)
        @test mr.get_time(snaps) == 2.3

        snaps = mr.evolve(snaps, simulation, 3.4)
        @test mr.get_time(snaps) == 3.4

        @test_throws ErrorException mr.evolve(snaps, simulation, 0.1)

        snaps = mr.evolve(mr.Spin(), simulation, 0.)
        @test mr.get_time(snaps) == 0.

        snaps = mr.evolve(mr.Spin(), simulation, 0., TR=2)
        @test mr.get_time(snaps) == 2.8

        snaps = mr.evolve(snaps, simulation, 0.5, TR=3)
        @test mr.get_time(snaps) == 6.1

        @test_throws ErrorException mr.evolve(snaps, simulation, 0.1)

        @test_throws MethodError mr.evolve(snaps, simulation)
    end
    @testset "Basic diffusion has no effect in constant fields" begin
        sequence = mr.read_pulseq(joinpath(@__DIR__, "pulseq", "gradient_echo_TE_2.seq"))
        no_diff = mr.Simulation([sequence], diffusivity=0., R2=0.3)
        with_diff = mr.Simulation([sequence], diffusivity=1., R2=0.3)
        spin_no_diff = mr.evolve(mr.Spin(), no_diff, 2.).spins[1]
        spin_with_diff = mr.evolve(mr.Spin(), with_diff, 2.).spins[1]
        @test spin_no_diff.position == SA[0, 0, 0]
        @test spin_with_diff.position != SA[0, 0, 0]
        @test mr.orientation.(spin_with_diff.orientations) == mr.orientation.(spin_no_diff.orientations)
        @test mr.transverse(spin_no_diff) ≈ exp(-0.6)
        @test abs(mr.longitudinal(spin_no_diff)) < Float64(1e-6)
    end
    @testset "Basic diffusion run within sphere" begin
        sequence = mr.read_pulseq(joinpath(@__DIR__, "pulseq", "gradient_echo_TE_20.seq"))
        sphere = mr.Spheres(radius=1.)
        Random.seed!(12)
        diff = mr.Simulation(sequence, diffusivity=2., geometry=sphere)
        snaps = mr.readout([mr.Spin(), mr.Spin()], diff, 0:5:mr.SequenceParts.repetition_time(sequence), return_snapshot=true)
        @test size(snaps) == (5, )
        for snap in snaps
            @test length(snap.spins) == 2
            for spin in snap.spins
                @test norm(spin.position) < 1.
                @test length(spin.orientations) == 1
            end
        end
    end
    @testset "Liminal extracellular encounters" begin
        geometry = mr.fix(mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1., permeability=Inf))],
            extracellular_fraction=0.2,
        ))
        internal = mr.Geometries.Internal
        Random.seed!(1234)
        collision = internal.detect_intersection(
            geometry,
            SA[0., 0., 0.],
            SA[10., 0., 0.],
        )
        @test collision !== nothing
        @test 0 < collision.distance < 1
        @test !collision.inside

        simulation = mr.Simulation([], geometry=geometry, diffusivity=0.1, timestep=0.2)
        sampling = mr.spin_sampling(geometry, mr.BoundingBox(2.), 1.)
        @test any(isempty, getfield.(sampling, :isinside))
        @test any(!isempty, getfield.(sampling, :isinside))
        snapshot = mr.Snapshot(sampling, 0.)
        evolved = mr.evolve(snapshot, simulation, 0.2)
        @test length(evolved.spins) == length(snapshot.spins)
    end
    @testset "Fallback bounding box for unsupported geometry" begin
        liminal = mr.fix(mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1.))],
            extracellular_fraction=0.2,
        ))
        simulation = mr.Simulation([], geometry=liminal, diffusivity=0.)
        snapshot = mr.Snapshot(1000, simulation)
        @test length(snapshot) == 1000
        @test all(
            all(position .>= -500.) && all(position .<= 500.)
            for position in mr.position.(snapshot)
        )
    end
    @testset "Liminal compartment density remains stable" begin
        f_extra = 0.2
        number_fraction = 1.
        radius = 1.
        geometry = mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=radius, permeability=Inf))],
            extracellular_fraction=f_extra,
        )
        simulation = mr.Simulation([], geometry=geometry, diffusivity=0.1, timestep=0.2)
        Random.seed!(4321)
        snapshot = mr.Snapshot(1000, simulation)

        function compartment_counts(snapshot)
            counts = zeros(Int, 2)
            for spin in snapshot.spins
                isempty(spin.isinside) ?
                    (counts[1] += 1) : (counts[first(spin.isinside)[1][1] + 1] += 1)
            end
            counts
        end

        compartment_history = Vector{Vector{Int}}()
        push!(compartment_history, compartment_counts(snapshot))
        for time in 1:6
            snapshot = mr.evolve(snapshot, simulation, time * 0.2)
            push!(compartment_history, compartment_counts(snapshot))
        end
        mean_counts = vec(mean(reduce(hcat, compartment_history), dims=2))
        expected = [
            f_extra * length(snapshot),
            (1 - f_extra) * number_fraction * length(snapshot),
        ]
        @test mean_counts ≈ expected rtol=0.2
    end
    @testset "Run simulation with multiple sequences at once" begin
        sequences = [
            build_sequence([mr.PulseEvent(flip_angle=0, phase=0.), 2., :readout, 1.]),
            build_sequence([mr.PulseEvent(flip_angle=90, phase=0.), 2., :readout, 1.]),
            build_sequence([mr.PulseEvent(flip_angle=90, phase=0.), 1., :readout, 1.])
        ]
        all_snaps = mr.Simulation(sequences, diffusivity=1., R2=1.)

        readouts = mr.readout(mr.Spin(), all_snaps)
        @test size(readouts) == (3,)

        # check relaxation
        @test mr.transverse(readouts[1]) ≈ 0. atol=1e-12
        @test mr.transverse(readouts[2]) ≈ exp(-2.)
        @test mr.transverse(readouts[3]) ≈ exp(-1.)
        @test mr.longitudinal(readouts[1]) ≈ 1.
        @test mr.longitudinal(readouts[2]) ≈ 0. atol=1e-12
        @test mr.longitudinal(readouts[3]) ≈ 0. atol=1e-12
    end

    @testset "Test readout identification" begin
        seq = build_sequence([
            mr.PulseEvent(flip_angle=0, phase=0.), 
            2., 
            :readout, 
            1.,
            :readout, 
            1.
        ])
        @test mr.get_readouts(seq, 0.) == (false, [
            mr.IndexedReadout(2., 0, 1),
            mr.IndexedReadout(3., 0, 2)
        ])
        @test mr.get_readouts(seq, 2.) == (false, [
            mr.IndexedReadout(2., 0, 1),
            mr.IndexedReadout(3., 0, 2)
        ])
        @test mr.get_readouts(seq, 2., readouts=[20., 1., 3., 2.]) == (false, [
            mr.IndexedReadout(2., 0, 4),
            mr.IndexedReadout(3., 0, 3),
            mr.IndexedReadout(20., 0, 1)
        ])
        @test mr.get_readouts(seq, 0., nTR=2) == (true, [
            mr.IndexedReadout(2., 1, 1),
            mr.IndexedReadout(3., 1, 2),
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ])
        @test mr.get_readouts(seq, 1., nTR=2) == (true, [
            mr.IndexedReadout(2., 1, 1),
            mr.IndexedReadout(3., 1, 2),
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ])
        @test mr.get_readouts(seq, 0., skip_TR=1) == (true, [
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ])
        @test mr.get_readouts(seq, 2.5, skip_TR=0) == (true, [
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ])
        @test mr.get_readouts(seq, 2., skip_TR=0) == (true, [
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ])

        @test_throws ErrorException mr.get_readouts(seq, 2., skip_TR=0, readouts=[2., 30])
        
        @test mr.get_readouts(seq, 2., skip_TR=0, readouts=[2., 4.00001]) == (true, [
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(8., 2, 2),
        ])
    end

end
