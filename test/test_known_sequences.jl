@testset "Validate expected output from known sequences" begin
    @testset "Diffusion MRI sequences" begin
        @testset "Relation between b-value, q-value and diffusion time" begin
            @test all(mr.derive_qval_time(80.) .≈ (0., 40.))
            @test all(mr.derive_qval_time(80., diffusion_time=10) .≈ (0., 10.))
            @test all(mr.derive_qval_time(80., bval=3) .≈ (sqrt(3/40), 40.))
            @test all(mr.derive_qval_time(80., bval=3, diffusion_time=1.) .≈ (sqrt(3), 1.))
            @test all(mr.derive_qval_time(80., bval=3, qval=2.) .≈ (2., 0.75))
            @test all(mr.derive_qval_time(80., diffusion_time=3, qval=2.) .≈ (2., 3.))
            @test all(mr.derive_qval_time(80., qval=2.) .≈ (2., 40.))
            @test all(mr.derive_qval_time(80., qval=1., diffusion_time=1, bval=1.) .≈ (1., 1.))

            @test_throws AssertionError mr.derive_qval_time(80., qval=1., diffusion_time=2., bval=0.)
        end
        @testset "Perfect PGSE with no diffusion" begin
            nspins = 300
            spins = mr.Snapshot(nspins)
            TE = 80.
            sequence = mr.perfect_dwi(bval=2., TE=TE)
            readout = mr.Simulation(spins, [sequence])
            append!(readout, TE * 1.01)
            @test length(readout.readout[1]) == 1
            snap = readout.readout[1][1]
            @test snap.time == TE
            @test mr.transverse(snap) ≈ nspins rtol=1e-2
        end
        @testset "Perfect PGSE with free diffusion" begin
            nspins = 1000
            spins = mr.Snapshot(nspins)
            TE = 80.
            sequence = mr.perfect_dwi(bval=0.3, TE=TE)
            readout = mr.Simulation(spins, [sequence], diffusivity=0.5)
            append!(readout, TE * 1.01)
            snap = readout.readout[1][1]
            @test snap.time == TE
            @test mr.transverse(snap) ≈ nspins * exp(-0.15) rtol=0.05
        end
    end
end
