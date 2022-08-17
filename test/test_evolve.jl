@testset "Evolve a single spin fully" begin
    @testset "Empty environment and sequence" begin
        simulation = mr.Simulation(mr.Spin(), [mr.Sequence(2.8)], mr.Microstructure(), store_every=0.5)
        append!(simulation, 2.8)
        snaps = mr.get_sequence.(simulation.regular, 1)
        time = 0.
        for snap in snaps
            @test snap.time == time
            time += 0.5
            @test mr.orientation(snap) == SA[0., 0., 1.]
            @test mr.longitudinal(snap) == 1.
            @test mr.transverse(snap) == 0.
        end
        @test length(snaps) == 6

        simulation = mr.Simulation([mr.Spin(), mr.Spin()], [mr.Sequence(2.8)], mr.Microstructure(), store_every=0.5)
        append!(simulation, 2.8)
        snaps = mr.get_sequence.(simulation.regular, 1)
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
        simulation = mr.Simulation(mr.Spin(), [mr.Sequence([mr.RFPulse(flip_angle=90)], 2.8)], mr.Microstructure(), store_every=0.5)
        append!(simulation, 2.8)
        snaps = mr.get_sequence.(simulation.regular, 1)
        @test mr.orientation(snaps[1]) ≈ SA[0., 0., 1.]
        for snap in snaps[2:end]
            @test mr.orientation(snap) ≈ SA[0., 1., 0.]
        end
        @test length(snaps) == 6
    end
    @testset "Ensure data is stored at requested time" begin
        simulation = mr.Simulation(mr.Spin(), [mr.Sequence(2.8)], mr.Microstructure(), store_every=0.5)
        append!(simulation, 2.8)
        @test length(simulation.regular) == 6
        @test simulation.latest[end].time == Float(2.8)
    end
    @testset "Basic diffusion has no effect in constant fields" begin
        sequence = mr.Sequence([mr.RFPulse(flip_angle=90)], 2.)
        no_diff = mr.Simulation(mr.Spin(), [sequence], mr.Microstructure(R2=mr.field(0.3)), store_every=0.5)
        with_diff = mr.Simulation(mr.Spin(), [sequence], mr.Microstructure(diffusivity=mr.field(1.), R2=mr.field(0.3)), store_every=0.5)
        append!(no_diff, sequence.TR)
        append!(with_diff, sequence.TR)
        spin_no_diff = mr.get_sequence(no_diff.latest[end].spins[1], 1)
        spin_with_diff = mr.get_sequence(with_diff.latest[end].spins[1], 1)
        @test spin_no_diff.position == SA[0, 0, 0]
        @test spin_with_diff.position != SA[0, 0, 0]
        @test spin_with_diff.orientations == spin_no_diff.orientations
        @test mr.transverse(spin_no_diff) ≈ exp(-0.6)
        @test abs(mr.longitudinal(spin_no_diff)) < Float(1e-6)
    end
    @testset "Basic diffusion changes spin orientation in spatially varying field" begin
        sequence = mr.Sequence([mr.RFPulse(flip_angle=90)], 2.)
        no_diff = mr.Simulation(mr.Spin(), [sequence], mr.Microstructure(R2=mr.field(0.3)), store_every=0.5)
        with_diff = mr.Simulation(mr.Spin(), [sequence], mr.Microstructure(diffusivity=mr.field(1.), R2=mr.field(SA[1., 0., 0.], 0.3)), store_every=0.5)
        with_diff_no_grad = mr.Simulation(mr.Spin(), [sequence], mr.Microstructure(diffusivity=mr.field(1.), R2=mr.field(0.3)), store_every=0.5)
        for res in (no_diff, with_diff, with_diff_no_grad)
            append!(res, sequence.TR)
        end
        spin_no_diff = mr.get_sequence(no_diff.latest[end].spins[1], 1)
        spin_with_diff = mr.get_sequence(with_diff.latest[end].spins[1], 1)
        spin_with_diff_no_grad = mr.get_sequence(with_diff_no_grad.latest[end].spins[1], 1)
        @test spin_no_diff.position == SA[0, 0, 0]
        @test spin_with_diff.position != SA[0, 0, 0]
        @test spin_with_diff_no_grad.position != SA[0, 0, 0]
        @test spin_with_diff.orientations != spin_no_diff.orientations
        @test spin_with_diff_no_grad.orientations == spin_no_diff.orientations
        @test abs(mr.longitudinal(spin_no_diff)) < 1e-6
    end
    @testset "Basic diffusion run within sphere" begin
        sequence = mr.Sequence([mr.RFPulse(flip_angle=90)], 2.)
        sphere = mr.Sphere(1.)
        Random.seed!(12)
        diff = mr.Simulation([mr.Spin(), mr.Spin()], [mr.Sequence(20.)], mr.Microstructure(diffusivity=mr.field(2.), geometry=sphere), store_every=0.5)
        append!(diff, sequence.TR)
        for snap in diff.regular
            @test length(snap.spins) == 2
            for spin in snap.spins
                @test norm(spin.position) < 1.
                @test length(spin.orientations) == 1
            end
        end
    end
    @testset "Run simulation with multiple sequences at once" begin
        sequences = [
            mr.Sequence([mr.RFPulse(flip_angle=0), mr.Readout(2.)], 3.),
            mr.Sequence([mr.RFPulse(flip_angle=90), mr.Readout(2.)], 3.),
            mr.Sequence([mr.RFPulse(flip_angle=90), mr.Readout(1.)], 2.),
        ]
        all_snaps = mr.Simulation(mr.Spin(), sequences, mr.Microstructure(diffusivity=mr.field(1.), R2=mr.field(1.)), store_every=0.2)
        append!(all_snaps, 3.)

        # check final time
        @test all_snaps.latest[end].time == 3.

        readouts = [r[1] for r in all_snaps.readout]

        # check relaxation
        @test mr.transverse(readouts[1]) ≈ 0. atol=1e-12
        @test mr.transverse(readouts[2]) ≈ exp(-2.)
        @test mr.transverse(readouts[3]) ≈ exp(-1.)
        @test mr.longitudinal(readouts[1]) ≈ 1.
        @test mr.longitudinal(readouts[2]) ≈ 0. atol=1e-12
        @test mr.longitudinal(readouts[3]) ≈ 0. atol=1e-12
    end
end
