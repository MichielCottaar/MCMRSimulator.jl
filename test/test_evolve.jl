@testset "test_evolve.jl" begin
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
            @test mr.orientation(snap) ≈ SA[0., -1., 0.]
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
        snaps = mr.readout([mr.Spin(), mr.Spin()], diff, 0:5:variables.duration(sequence), return_snapshot=true)
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
        seq = build_sequence() do 
            Sequence([
                InstantPulse(flip_angle=0, phase=0.), 
                2., 
                SingleReadout(), 
                1.,
                SingleReadout(), 
                1.
            ]) 
        end
        @test collect(mr.get_readouts(seq, 0.)) == [
            mr.IndexedReadout(2., 1, 1),
            mr.IndexedReadout(3., 1, 2)
        ]
        @test collect(mr.get_readouts(seq, 2.)) == [
            mr.IndexedReadout(3., 1, 2)
        ]
        @test collect(mr.get_readouts(seq, 2., readouts=[20., 1., 3., 2.])) == [
            mr.IndexedReadout(3., 1, 3)
            mr.IndexedReadout(20., 1, 1)
        ]
        @test collect(mr.get_readouts(seq, 0., nTR=2)) == [
            mr.IndexedReadout(2., 1, 1),
            mr.IndexedReadout(3., 1, 2),
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ]
        @test collect(mr.get_readouts(seq, 1., nTR=2)) == [
            mr.IndexedReadout(2., 1, 1),
            mr.IndexedReadout(3., 1, 2),
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ]
        @test collect(mr.get_readouts(seq, 0., skip_TR=1)) == [
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ]
        @test collect(mr.get_readouts(seq, 2.5, skip_TR=0)) == [
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ]
        @test collect(mr.get_readouts(seq, 2., skip_TR=0)) == [
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(7., 2, 2),
        ]

        @test_throws ErrorException mr.get_readouts(seq, 2., skip_TR=0, readouts=[2., 30])
        
        @test collect(mr.get_readouts(seq, 2., skip_TR=0, readouts=[2., 4.00001])) == [
            mr.IndexedReadout(6., 2, 1),
            mr.IndexedReadout(8., 2, 2),
        ]
    end

end
