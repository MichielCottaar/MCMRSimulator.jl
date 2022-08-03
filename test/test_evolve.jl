@testset "Evolve a single spin fully" begin
    @testset "Empty environment and sequence" begin
        snaps = evolve_TR(Spin(), Sequence(2.8), Microstructure(), store_every=0.5).regular
        time = 0.
        for snap in snaps
            @test snap.time == time
            time += 0.5
            @test vector(snap) == SA_F64[0., 0., 1.]
            @test longitudinal(snap) == 1.
            @test transverse(snap) == 0.
        end
        @test length(snaps) == 5

        snaps = evolve_TR([Spin(), Spin()], Sequence(2.8), Microstructure(), store_every=0.5).regular
        time = 0.
        for snap in snaps
            @test snap.time == time
            time += 0.5
            @test vector(snap) == SA_F64[0., 0., 2.]
            @test longitudinal(snap) == 2.
            @test transverse(snap) == 0.
        end
        @test length(snaps) == 5
    end
    @testset "Gradient echo sequence" begin
        snaps = evolve_TR(Spin(), Sequence([RFPulse(flip_angle=90)], 2.8), Microstructure(), store_every=0.5).regular
        for snap in snaps
            @test vector(snap) ≈ SA_F64[0., 1., 0.]
        end
        @test length(snaps) == 5
    end
    @testset "Ensure data is stored at final TR" begin
        snaps = evolve_TR(Spin(), Sequence(2.), Microstructure(), store_every=0.5)
        @test length(snaps) == 5
    end
    @testset "Basic diffusion has no effect in constant fields" begin
        sequence = Sequence([RFPulse(flip_angle=90)], 2.)
        no_diff = evolve_TR(Spin(), sequence, Microstructure(R2=field(0.3)), store_every=0.5)
        with_diff = evolve_TR(Spin(), sequence, Microstructure(diffusivity=field(1.), R2=field(0.3)), store_every=0.5)
        with_diff_grad = evolve_TR(Spin(), sequence, Microstructure(diffusivity=field(1.), R2=field(0.3)), store_every=0.5)
        spin_no_diff = no_diff[end].spins[1]
        spin_with_diff = with_diff[end].spins[1]
        @test spin_no_diff.position == SA_F64[0, 0, 0]
        @test spin_with_diff.position != SA_F64[0, 0, 0]
        @test spin_with_diff.orientation == spin_no_diff.orientation
    end
    @testset "Basic diffusion changes spin orientation in spatially varying field" begin
        sequence = Sequence([RFPulse(flip_angle=90)], 2.)
        no_diff = evolve_TR(Spin(), sequence, Microstructure(R2=field(SA_F64[1., 0, 0], 0.3)), store_every=0.5)
        with_diff = evolve_TR(Spin(), sequence, Microstructure(diffusivity=field(1.), R2=field(SA_F64[1., 0., 0.], 0.3)), store_every=0.5)
        with_diff_no_grad = evolve_TR(Spin(), sequence, Microstructure(diffusivity=field(1.), R2=field(0.3)), store_every=0.5)
        spin_no_diff = no_diff[end].spins[1]
        spin_with_diff = with_diff[end].spins[1]
        spin_with_diff_no_grad = with_diff_no_grad[end].spins[1]
        @test spin_no_diff.position == SA_F64[0, 0, 0]
        @test spin_with_diff.position != SA_F64[0, 0, 0]
        @test spin_with_diff_no_grad.position != SA_F64[0, 0, 0]
        @test spin_with_diff.orientation != spin_no_diff.orientation
        @test spin_with_diff_no_grad.orientation == spin_no_diff.orientation
    end
    @testset "Basic diffusion run within sphere" begin
        sequence = Sequence([RFPulse(flip_angle=90)], 2.)
        sphere = Sphere(1.)
        Random.seed!(12)
        diff = evolve_TR([Spin(), Spin()], Sequence(20.), Microstructure(diffusivity=field(2.), geometry=sphere), store_every=0.5)
        for snap in diff
            @test length(snap) == 2
            for spin in snap
                @test norm(spin.position) < 1.
            end
        end
    end
    @testset "Run simulation with multiple sequences at once" begin
        sequences = [
            Sequence([RFPulse(flip_angle=0)], 2.),
            Sequence([RFPulse(flip_angle=90)], 2.),
            Sequence([RFPulse(flip_angle=90)], 1.),
        ]
        all_snaps = evolve_TR(Spin(), sequences, Microstructure(diffusivity=field(1.), R2=field(1.)), store_every=0.2)

        # check final time
        @test all_snaps[1][end].time == 2.
        @test all_snaps[2][end].time == 2.
        @test all_snaps[3][end].time == 1.

        # check relaxation
        @test transverse(all_snaps[1][end]) ≈ 0. atol=1e-12
        @test transverse(all_snaps[2][end]) ≈ exp(-2.)
        @test transverse(all_snaps[3][end]) ≈ exp(-1.)
        @test longitudinal(all_snaps[1][end]) ≈ 1.
        @test longitudinal(all_snaps[2][end]) ≈ 0. atol=1e-12
        @test longitudinal(all_snaps[3][end]) ≈ 0. atol=1e-12

        # Compare spin positions after a little bit of diffusion
        ref = all_snaps[1][4]
        @test time(ref) > 0
        @test ref[1].position == all_snaps[2][4][1].position
        @test ref[1].position == all_snaps[3][4][1].position
    end
end
