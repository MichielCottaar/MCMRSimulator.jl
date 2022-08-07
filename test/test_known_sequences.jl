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
            spins = [mr.Spin(position=randn(mr.PosVector) * 100.) for _ in 1:100]
            sequence = mr.perfect_dwi(bval=2., TE=80.)
            micro = mr.Microstructure()
            readout = mr.Simulation(spins, [sequence], micro)
            push!(readout, 90.)
            @test length(readout.readout[sequence][1]) == 1
            snap = readaout.readaout[sequence][1]
            @test snap.time == 80.
            @test transverse(snap) ≈ 100. rtol=1e-2
        end
        if false
            @testset "Perfect PGSE with free diffusion" begin
                spins = [mr.Spin(position=randn(mr.PosVector) * 100.) for _ in 1:100]
                sequence = mr.perfect_dwi(bval=2.)
                micro = mr.Microstructure(diffusivity=mr.field(3.))
                readout = mr.evolve_TR(spins, sequence, micro)
                println(log(transverse(readout.data[1])))
                @test transverse(readout.data[1]) ≈ 100. * exp(-6.) rtol=0.1
            end
        end
    end
end
