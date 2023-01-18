@testset "Evolve a single spin fully" begin
    @testset "Timings of readouts" begin
        @testset "Explicitly setting timestep" begin
            s1 = mr.dwi(bval=0.)
            s2 = mr.dwi(bval=1.)

            @test all(mr.get_times(mr.Simulation([s1, s2], timestep=0.5), 0., 80.) .== 0:0.5:80.)
            @test all(mr.get_times(mr.Simulation([s1, s2], timestep=0.501), 0., 80.) .== 0:0.5:80.)
            @test all(mr.get_times(mr.Simulation([s1, s2], timestep=0.501), 0., 100.) .== 0:0.5:100.)
        end
        @testset "Setting gradient_precision" begin
            s1 = mr.dwi(bval=0.)
            s2 = mr.dwi(bval=1.)

            kwargs = Dict(
                :gradient_precision => 0.01,
                :sample_displacement => 0,
                :sample_off_resonance => 0,
            )

            @test all(mr.get_times(mr.Simulation(s1; kwargs...), 0., 80.) .== [0., 40., 80.])
            @test all(mr.get_times(mr.Simulation(s1; kwargs...), 90., 120.) .== [90., 120.])
            @test all(mr.get_times(mr.Simulation([s1, s2]; kwargs...), 0., 80.) .== [0., 40., 80.])
            @test all(mr.get_times(mr.Simulation([s1, s2]; kwargs...), 90., 120.) .== [90., 120.])

            kwargs[:diffusivity] = 1.
            @test all(mr.get_times(mr.Simulation(s1; kwargs...), 0., 80.) .== [0., 40., 80.])
            @test all(mr.get_times(mr.Simulation(s1; kwargs...), 90., 120.) .== [90., 120.])
            @test length(mr.get_times(mr.Simulation([s1, s2]; kwargs...), 0., 80.)) > 3
            @test length(union(mr.get_times(mr.Simulation([s1, s2]; kwargs...), 0., 80.), [0., 40., 80.])) == 3.
            @test all(mr.get_times(mr.Simulation([s1, s2]; kwargs...), 90., 120.) .== [90., 120.])
        end
        @testset "Setting sample_displacement" begin
            geometry = mr.Annulus(1., 2., myelin=true)
            s1 = mr.dwi(bval=0.)
            s2 = mr.dwi(bval=1., readout_time=40)

            kwargs = Dict(
                :gradient_precision => 0,
                :sample_displacement => 4,
                :sample_off_resonance => 0,
                :geometry => geometry,
                :diffusivity => 1.,
            )

            @test all(mr.get_times(mr.Simulation(s1; kwargs...), 0., 80.) .== 0:10:80)
            @test all(mr.get_times(mr.Simulation(s1; kwargs...), 90., 130.) .== 90:10:130)

            @test all(mr.get_times(mr.Simulation([s1, s2]; kwargs...), 0., 80.) .== [0:5:60..., 70, 80])
            @test all(mr.get_times(mr.Simulation([s1, s2]; kwargs...), 90., 130.) .== 90:10:130)
        end
        @testset "Setting sample_off_resonance" begin
            geometry = mr.Annulus(1., 2., myelin=true)
            s1 = mr.dwi(bval=0.)
            s2 = mr.dwi(bval=1., readout_time=40)

            kwargs = Dict(
                :gradient_precision => 0,
                :sample_displacement => 0,
                :sample_off_resonance => 4,
                :geometry => geometry,
                :diffusivity => 1.,
            )

            for sequences in (s1, [s1, s2])
                @test all(mr.get_times(mr.Simulation(sequences; kwargs...), 0., 80.) .== 0:10:80)
                @test all(mr.get_times(mr.Simulation(sequences; kwargs...), 85., 125.) .== 85:20:125)
            end
        end
    end
    @testset "Empty environment and sequence" begin
        simulation = mr.Simulation(mr.Sequence(TR=2.8))
        snaps = mr.trajectory(zeros(3), simulation, 0:0.5:2.8)
        time = 0.
        for snap in snaps
            @test snap.time == time
            time += 0.5
            @test mr.orientation(snap) == SA[0., 0., 1.]
            @test mr.longitudinal(snap) == 1.
            @test mr.transverse(snap) == 0.
        end
        @test length(snaps) == 6

        simulation = mr.Simulation(mr.Sequence(TR=2.8))
        snaps= mr.trajectory([mr.Spin(), mr.Spin()], simulation, 0:0.5:2.8)
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
        simulation = mr.Simulation(mr.Sequence(pulses=[mr.RFPulse(flip_angle=90)], TR=2.8))
        snaps = mr.trajectory(zeros(3), simulation, 0:0.5:2.8)
        @test mr.orientation(snaps[1]) ≈ SA[0., 0., 1.]
        for snap in snaps[2:end]
            @test mr.orientation(snap) ≈ SA[0., 1., 0.]
        end
        @test length(snaps) == 6
    end
    @testset "Ensure data is stored at requested time" begin
        simulation = mr.Simulation([mr.Sequence(TR=2.8)])
        snaps = mr.trajectory(mr.Spin(), simulation, 0:0.5:2.8)
        @test length(snaps) == 6

        snaps = mr.evolve(mr.Spin(), simulation, 2.3)
        @test time(snaps) == 2.3

        snaps = mr.evolve(snaps, simulation)
        @test time(snaps) == 2.8

        snaps = mr.evolve(mr.Spin(), simulation)
        @test time(snaps) == 2.8
    end
    @testset "Basic diffusion has no effect in constant fields" begin
        sequence = mr.Sequence(pulses=[mr.RFPulse(flip_angle=90)], TR=2.)
        no_diff = mr.Simulation([sequence], R2=mr.field(0.3))
        with_diff = mr.Simulation([sequence], diffusivity=1., R2=mr.field(0.3))
        spin_no_diff = mr.evolve(mr.Spin(), no_diff).spins[1]
        spin_with_diff = mr.evolve(mr.Spin(), with_diff).spins[1]
        @test spin_no_diff.position == SA[0, 0, 0]
        @test spin_with_diff.position != SA[0, 0, 0]
        @test spin_with_diff.orientations == spin_no_diff.orientations
        @test mr.transverse(spin_no_diff) ≈ exp(-0.6)
        @test abs(mr.longitudinal(spin_no_diff)) < Float(1e-6)
    end
    @testset "Basic diffusion changes spin orientation in spatially varying field" begin
        sequence = mr.Sequence(pulses=[mr.RFPulse(flip_angle=90)], TR=2.)
        no_diff = mr.Simulation([sequence], mr.Microstructure(R2=mr.field(0.3)))
        with_diff = mr.Simulation([sequence], mr.Microstructure(diffusivity=mr.field(1.), R2=mr.field([1., 0., 0.], 0.3)))
        with_diff_no_grad = mr.Simulation([sequence], mr.Microstructure(diffusivity=mr.field(1.), R2=mr.field(0.3)))

        spin_no_diff = mr.evolve(mr.Spin(), no_diff).spins[1]
        spin_with_diff = mr.evolve(mr.Spin(), with_diff).spins[1]
        spin_with_diff_no_grad = mr.evolve(mr.Spin(), with_diff_no_grad).spins[1]

        @test spin_no_diff.position == SA[0, 0, 0]
        @test spin_with_diff.position != SA[0, 0, 0]
        @test spin_with_diff_no_grad.position != SA[0, 0, 0]
        @test spin_with_diff.orientations != spin_no_diff.orientations
        @test spin_with_diff_no_grad.orientations == spin_no_diff.orientations
        @test abs(mr.longitudinal(spin_no_diff)) < 1e-6
    end
    @testset "Basic diffusion run within sphere" begin
        sequence = mr.Sequence(pulses=[mr.RFPulse(flip_angle=90)], TR=2.)
        sphere = mr.Sphere(1.)
        Random.seed!(12)
        diff = mr.Simulation([mr.Sequence(TR=20.)], diffusivity=2., geometry=sphere)
        snaps = mr.trajectory([mr.Spin(), mr.Spin()], diff, 0:0.5:sequence.TR)
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
            mr.Sequence(pulses=[mr.RFPulse(flip_angle=0), mr.Readout(2.)], TR=3.),
            mr.Sequence(pulses=[mr.RFPulse(flip_angle=90), mr.Readout(2.)], TR=3.),
            mr.Sequence(pulses=[mr.RFPulse(flip_angle=90), mr.Readout(1.)], TR=2.),
        ]
        all_snaps = mr.Simulation(sequences, diffusivity=1., R2=mr.field(1.))

        readouts = [r[1] for r in mr.readout(mr.Spin(), all_snaps)]

        # check relaxation
        @test mr.transverse(readouts[1]) ≈ 0. atol=1e-12
        @test mr.transverse(readouts[2]) ≈ exp(-2.)
        @test mr.transverse(readouts[3]) ≈ exp(-1.)
        @test mr.longitudinal(readouts[1]) ≈ 1.
        @test mr.longitudinal(readouts[2]) ≈ 0. atol=1e-12
        @test mr.longitudinal(readouts[3]) ≈ 0. atol=1e-12
    end
end
