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
            sequence = mr.perfect_dwi(bval=2., TE=80.)
            sim = mr.Simulation([sequence])
            readout = mr.readout(nspins, sim)
            @test length(readout[1]) == 1
            snap = readout[1][1]
            @test snap.time == 80.
            @test mr.transverse(snap) ≈ nspins rtol=1e-2
        end
        @testset "Perfect PGSE with free diffusion" begin
            nspins = 1000
            TE = 80.
            sequence = mr.perfect_dwi(bval=0.3, TE=TE)
            sim = mr.Simulation([sequence], diffusivity=0.5)
            readout = mr.readout(nspins, sim)
            snap = readout[1][1]
            @test snap.time == TE
            @test mr.transverse(snap) ≈ nspins * exp(-0.15) rtol=0.05
        end
        @testset "Diffusion within the sphere" begin
            for radius in (0.5, 1.)
                sphere = mr.Sphere(radius)
                all_pos = [rand(3) .* 2 .- 1 for _ in 1:30000]
                snap = mr.Snapshot([mr.Spin(position=pos .* radius) for pos in all_pos if norm(pos) < 1.])
                @testset "Stejskal-Tanner approximation at long diffusion times" begin
                    # equation 4 from Balinov, B. et al. (1993) ‘The NMR Self-Diffusion Method Applied to Restricted Diffusion. Simulation of Echo Attenuation from Molecules in Spheres and between Planes’, Journal of Magnetic Resonance, Series A, 104(1), pp. 17–25. doi:10.1006/jmra.1993.1184.
                    qvals = [0.01, 0.1, 1.]
                    sequences = [mr.perfect_dwi(TE=101, diffusion_time=100, qval=qval) for qval in qvals]
                    simulation = mr.Simulation(sequences; geometry=sphere, diffusivity=3.)
                    at_readout = mr.readout(snap, simulation)
                    for (qval, readouts) in zip(qvals, at_readout)
                        readout = readouts[1]
                        factor = qval * radius
                        expected = 9 * (factor * cos(factor) - sin(factor)) ^2 / factor^6
                        @test readout.time == 101
                        @test mr.transverse(readout) ≈ length(readout.spins) * expected rtol=0.05
                    end
                end
            end
        end
        @testset "Diffusion between two planes" begin
            for distance in [0.5, 1]
                walls = mr.walls(positions=[0., distance])
                snap = mr.Snapshot([mr.Spin(position=rand(3) * distance) for _ in 1:2000])
                @testset "Stejskal-Tanner approximation at long diffusion times for a=$distance" begin
                    # equation 6 from Balinov, B. et al. (1993) ‘The NMR Self-Diffusion Method Applied to Restricted Diffusion. Simulation of Echo Attenuation from Molecules in Spheres and between Planes’, Journal of Magnetic Resonance, Series A, 104(1), pp. 17–25. doi:10.1006/jmra.1993.1184.
                    qvals = [0.01, 0.1, 1.]
                    sequences = [mr.perfect_dwi(TE=101, diffusion_time=100, qval=qval, orientation=[1., 0., 0.]) for qval in qvals]
                    simulation = mr.Simulation(sequences; geometry=walls, diffusivity=3.)
                    at_readout = mr.readout(snap, simulation)
                    for (qval, readouts) in zip(qvals, at_readout)
                        readout = readouts[1]
                        factor = qval * distance
                        #factor = 2 * π * qval * distance
                        expected = 2 * (1 - cos(factor)) / factor^2
                        @test readout.time == 101
                        @test mr.transverse(readout) ≈ length(readout.spins) * expected rtol=0.05
                    end
                end
            end
        end
    end
end
