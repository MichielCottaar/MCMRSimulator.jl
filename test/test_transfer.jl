@testset "Test magnetisation transfer" begin
    @testset "Test signal loss between two walls" begin
        "Fraction of axons hitting the wall"
        frachit(x) = ((1 - erf(x)) + (1 - exp(-x*x)) / (sqrt(π) * x)) / 2
        frachit(wall_dist, diffusivity, timestep) = frachit(wall_dist / (2 * sqrt(diffusivity * timestep)))

        function test_MT_walls(wall_dist, diffusivity, timestep; nspins=100000, transfer=0.5)
            Random.seed!(1234)
            geometry = mr.walls(MT_fraction=transfer)
            spins = [mr.Spin(position=Random.rand(3) .* wall_dist) for _ in 1:nspins]
            sequence = mr.Sequence(pulses=[mr.RFPulse(flip_angle=90)], TR=1e5)
            simulation = mr.Simulation(sequence, geometry=geometry, diffusivity=diffusivity, timestep=timestep)
            signal = mr.signal(spins, simulation, timestep)
            @test length(signal) == 2
            fhit = frachit(wall_dist, diffusivity, timestep)
            @test transfer * fhit ≈ (1 - mr.transverse(signal[end]) / nspins) rtol=0.03
        end
        test_MT_walls(2., 3., 0.1; transfer=0.5)
        test_MT_walls(2., 3., 0.1; transfer=0.1)
        test_MT_walls(10., 3., 0.1; transfer=0.5)
        test_MT_walls(10., 1., 0.1; transfer=0.5)
    end
end
