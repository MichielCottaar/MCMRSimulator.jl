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
    @testset "Basic diffusion run within sphere" begin
        sequence = Sequence([RFPulse(flip_angle=90)], 2.)
        sphere = Sphere(1.)
        diff = evolve([Spin(), Spin()], Microstructure(diffusivity=field(2.), geometry=sphere), Sequence(20.), yield_every=0.5)
        for snap in diff
            @test length(snap) == 2
            for spin in snap
                @test norm(spin.position) < 1.
            end
        end
    end
end
