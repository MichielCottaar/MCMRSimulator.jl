@testset "test_known_sequences.jl: Validate expected output from known sequences" begin
    @testset "Diffusion MRI sequences" begin
        @testset "PGSE with various gradient durations and no diffusion" begin
            for (δ, Δ) in [
                (nothing, 30.),
                (10., 30.),
                (nothing, nothing),
                (10., nothing),
                (0, nothing),
                (0, 70.),
                (0, 0.1),
            ]
                @testset "δ=$δ; Δ=$Δ" begin
                    nspins = 300
                    gradient = (isnothing(δ) || !iszero(δ)) ? (type=:trapezoid, δ=δ) : (type=:instant, )
                    sequence = @test_nowarn DWI(bval=0.3, TE=80., Δ=Δ, gradient=gradient, scanner=Siemens_Connectom)
                    sim = mr.Simulation(sequence, diffusivity=0.)
                    snap = mr.readout(nspins, sim, return_snapshot=true)
                    @test snap.time ≈ 80.
                    @test mr.transverse(snap) ≈ nspins rtol=1e-2
                end
            end
        end
        @testset "PGSE with various gradient durations and free diffusion" begin
            for (δ, Δ) in [
                (nothing, 30.),
                (10, 30.),
                (20, 40.),
                (nothing, nothing),
                (10, nothing),
                (0, nothing),
                (0, 70),
                (0, 30.),
            ]
                nspins = 3000
                TE = 80.
                gradient = (isnothing(δ) || !iszero(δ)) ? (type=:trapezoid, δ=δ) : (type=:instant, )
                sequence = DWI(bval=0.3, TE=TE, Δ=Δ, gradient=gradient)
                sim = mr.Simulation(sequence, diffusivity=0.5)
                snap = mr.readout(nspins, sim, return_snapshot=true)
                @test snap.time ≈ TE
                @test mr.transverse(snap) ≈ nspins * exp(-0.15) rtol=0.05
            end
        end
        @testset "PGSE in realistic scanner with free diffusion" begin
            nspins = 30000
            sequence = DWI(bval=2, TE=80., scanner=Siemens_Prisma)
            sim = mr.Simulation(sequence, diffusivity=0.5)
            snap = mr.readout(nspins, sim)
            @test mr.transverse(snap) ≈ nspins * exp(-1) rtol=0.05
        end
        @testset "Diffusion within the sphere" begin
            for radius in (0.5, 1.)
                sphere = mr.Spheres(radius=radius)
                all_pos = [rand(3) .* 2 .- 1 for _ in 1:30000]
                snap = mr.Snapshot([mr.Spin(position=pos .* radius) for pos in all_pos if norm(pos) < 1.])
                @testset "Stejskal-Tanner approximation at long diffusion times" begin
                    # equation 4 from Balinov, B. et al. (1993) ‘The NMR Self-Diffusion Method Applied to Restricted Diffusion. Simulation of Echo Attenuation from Molecules in Spheres and between Planes’, Journal of Magnetic Resonance, Series A, 104(1), pp. 17–25. doi:10.1006/jmra.1993.1184.
                    qvals = [0.01, 0.1, 1.]
                    sequences = [DWI(TE=101, diffusion_time=100, qval=qval, gradient=(type=:instant, )) for qval in qvals]
                    simulation = mr.Simulation(sequences; geometry=sphere, diffusivity=3.)
                    at_readout = mr.readout(snap, simulation, return_snapshot=true)
                    for (qval, readout) in zip(qvals, at_readout)
                        factor = qval * radius
                        expected = 9 * (factor * cos(factor) - sin(factor)) ^2 / factor^6
                        @test readout.time == 101
                        @test log(mr.transverse(readout) / length(readout.spins)) ≈ log(expected) rtol=0.05
                    end
                end
            end
        end
        @testset "Diffusion between two planes" begin
            for distance in [0.5, 1]
                walls = mr.Walls(position=[0., distance])
                Random.seed!(1234)
                snap = mr.Snapshot([mr.Spin(position=rand(3) * distance) for _ in 1:5000])
                @testset "Stejskal-Tanner approximation at long diffusion times for a=$distance" begin
                    # equation 6 from Balinov, B. et al. (1993) ‘The NMR Self-Diffusion Method Applied to Restricted Diffusion. Simulation of Echo Attenuation from Molecules in Spheres and between Planes’, Journal of Magnetic Resonance, Series A, 104(1), pp. 17–25. doi:10.1006/jmra.1993.1184.
                    qvals = [0.01, 0.1, 1.]
                    sequences = [DWI(TE=101, diffusion_time=100, qval=qval, gradient=(type=:instant, orientation=[1., 0., 0.])) for qval in qvals]
                    simulation = mr.Simulation(sequences; geometry=walls, diffusivity=3.)
                    at_readout = mr.readout(snap, simulation, return_snapshot=true)
                    for (qval, readout) in zip(qvals, at_readout)
                        factor = qval * distance
                        #factor = 2 * π * qval * distance
                        expected = 2 * (1 - cos(factor)) / factor^2
                        @test readout.time == 101
                        @test log(mr.transverse(readout) / length(snap)) ≈ log(expected) rtol=0.05
                    end
                end
                @testset "Mitra approximation at long diffusion times" begin
                    # equation 3 from Mitra, P.P. et al. (1992) ‘Diffusion propagator as a probe of the structure of porous media’, Physical Review Letters, 68(24), pp. 3555–3558. doi:10.1103/physrevlett.68.3555.
                    diffusion_times = [0.003, 0.01]
                    sequences = [
                        DWI(diffusion_time=dt, TE=2, bval=2.)
                        for dt in diffusion_times
                    ]
                    simulation = mr.Simulation(sequences; geometry=walls, diffusivity=1.)
                    at_readout = mr.readout(snap, simulation)

                    for (dt, readout_dt) in zip(diffusion_times, at_readout)
                        effective_diffusion = 1. - 4 / 3 * sqrt(π * dt) / distance
                        signal = mr.transverse(readout_dt)
                        @test log(signal / length(snap)) ≈ -2. * effective_diffusion rtol=0.2
                    end
                end
            end
        end
    end
end
