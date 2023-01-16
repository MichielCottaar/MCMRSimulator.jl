@testset "Test magnetisation transfer" begin
    @testset "Test signal loss between two walls" begin
        "Fraction of axons hitting the wall"
        frachit(x) = ((1 - erf(x)) + (1 - exp(-x*x)) / (sqrt(π) * x)) / 2
        frachit(wall_dist, diffusivity, timestep) = frachit(wall_dist / (2 * sqrt(diffusivity * timestep)))

        function test_MT_walls(wall_dist, diffusivity, timestep; nspins=100000, transfer=0.5)
            Random.seed!(1234)
            actual_transfer = 1 - (1 - transfer) ^ (1/sqrt(timestep))
            geometry = mr.walls(MT_fraction=actual_transfer)
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
    @testset "Test that transfer rate does not depend on timestep" begin
        Random.seed!(1234)
        geometry = mr.walls(repeats=1, MT_fraction=0.1)
        sequence = mr.Sequence(pulses=[mr.RFPulse(flip_angle=90)], TR=1e5)

        reference = nothing
        for timestep in (0.01, 0.1, 1)
            simulation = mr.Simulation(sequence, geometry=geometry, diffusivity=1., timestep=timestep)
            signal = mean([mr.transverse(mr.evolve(10000, simulation, 10.)) for _ in 1:Int(timestep/0.01)]) / 10000
            @show (timestep, log(signal))
            if isnothing(reference)
                reference = signal
            else
                @test log(reference) ≈ log(signal) rtol=0.1
            end
        end
    end
end
