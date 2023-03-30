@testset "test_transfer.jl" begin
@testset "Generate stuck particles" begin
    nspins = Int(1e5)
    @testset "Particles on repeating wall" begin
        s = mr.Simulation([], geometry=mr.walls(repeats=1., rotation=[1, 1, 0]), surface_density=1, dwell_time=1, diffusivity=1)
        all_spins = mr.Snapshot(nspins, s, 5.)
        @test !any([mr.isinside(s.geometry, p) > 0 for p in all_spins])
        free_spins = filter(s->!mr.stuck(s), [all_spins...])
        for dim in 1:3
            arr = [mr.position(s)[dim] for s in free_spins]
            @test mean(arr) ≈ 0 atol=0.1
            @test maximum(abs.(arr)) < 5
            @test maximum(abs.(arr)) > 4.9
        end

        stuck_spins = filter(mr.stuck, [all_spins...])
        @test length(stuck_spins) ≈ nspins / 2 rtol=0.1
        (x, y, z) = [[mr.position(s)[dim] for s in stuck_spins] for dim in 1:3]
        for arr in (x, y, z)
            @test mean(arr) ≈ 0 atol=0.1
            @test maximum(abs.(arr)) < 5
            @test maximum(abs.(arr)) > 4.9
        end
        limited_dim = (x + y) / sqrt(2)
        @test maximum(limited_dim) ≈ 7
        @test minimum(limited_dim) ≈ -7
    end
    @testset "Particles on two annuli" begin
        s = mr.Simulation([], geometry=mr.annuli(0.4, 0.6, positions=[[0, 0], [1, 1]], surface_density=[1, 2], dwell_time=1), diffusivity=1)
        all_spins = mr.Snapshot(nspins, s, 2)
        stuck_spins = filter(mr.stuck, all_spins)

        relative_to_free = 2π * (0.4 + 0.6) * (1 + 2) / 4^2  # surface area * surface_density / volume
        fraction_stuck = relative_to_free / (1 + relative_to_free)
        @test fraction_stuck ≈ length(stuck_spins) / nspins rtol=0.1

        (x, y, z) = [[mr.position(s)[dim] for s in stuck_spins] for dim in 1:3]

        # analyse annulus at origin
        radius = @. sqrt(x^2 + y^2)
        inner_cylinder = radius .< 0.5
        @test sum(inner_cylinder) / nspins ≈ fraction_stuck / 3 * 0.4 rtol=0.1
        @test all(radius[inner_cylinder] .≈ 0.4)
        @test mean([mr.isinside(s.geometry, p) for p in stuck_spins[inner_cylinder]]) ≈ 1.5 rtol=0.1
        outer_cylinder = 0.5 .< radius .< 0.7
        @test sum(outer_cylinder) / nspins ≈ fraction_stuck / 3 * 0.6 rtol=0.1
        @test all(radius[outer_cylinder] .≈ 0.6)
        @test mean([mr.isinside(s.geometry, p) for p in stuck_spins[outer_cylinder]]) ≈ 0.5 rtol=0.1

        # same for shifted annulus
        radius = @. sqrt((x - 1)^2 + (y - 1)^2)
        inner_cylinder = radius .< 0.5
        @test sum(inner_cylinder) / nspins ≈ fraction_stuck * 2 / 3 * 0.4 rtol=0.1
        @test all(radius[inner_cylinder] .≈ 0.4)
        @test mean([mr.isinside(s.geometry, p) for p in stuck_spins[inner_cylinder]]) ≈ 1.5 rtol=0.1
        outer_cylinder = 0.5 .< radius .< 0.7
        @test sum(outer_cylinder) / nspins ≈ fraction_stuck * 2 / 3 * 0.6 rtol=0.1
        @test all(radius[outer_cylinder] .≈ 0.6)
        @test mean([mr.isinside(s.geometry, p) for p in stuck_spins[outer_cylinder]]) ≈ 0.5 rtol=0.1
    end
    @testset "Correct number of stuck particles for long timesteps" begin
        Random.seed!(1234)
        geometry = mr.walls(repeats=1, surface_density=0.5, dwell_time=1.)
        sim = mr.Simulation([], diffusivity=3., geometry=geometry, max_timestep=10)
        nspins = 10000
        snapshot = mr.evolve(nspins, sim, 0)
        fraction_stuck = sum(mr.stuck.(snapshot)) / nspins
        @test fraction_stuck ≈ 1//3 rtol=0.05
        snapshot2 = mr.evolve(snapshot, sim, 300)
        fraction_stuck = sum(mr.stuck.(snapshot2)) / nspins
        # TODO: this test is far less strict than it should be
        @test fraction_stuck < 0.5
        displacement = mr.position.(snapshot2) .- mr.position.(snapshot)
        for dim in (2, 3)
            @test var([d[dim] for d in displacement]) / (2 * 300) ≈ 2//3 * 3 rtol=0.05  # spends 2/3rd of time as free spin
        end
    end
    @testset "Particles on two repeating cylinders" begin
        s = mr.Simulation([], geometry=mr.cylinders(0.6, positions=[[0, 0], [1, 1]], repeats=[2, 2], surface_density=[1, 2]), diffusivity=1, dwell_time=1)
    
        all_spins = mr.Snapshot(nspins, s, 1)
        stuck_spins = filter(mr.stuck, all_spins)

        relative_to_free = 2π * 0.6 * (1 + 2) / 2^2  # surface area * surface_density / volume
        fraction_stuck = relative_to_free / (1 + relative_to_free)
        @test fraction_stuck ≈ length(stuck_spins) / nspins rtol=0.1

        (x, y, z) = [[mr.position(s)[dim] for s in stuck_spins] for dim in 1:3]
        radius = @. sqrt(x^2 + y^2)
        origin_cylinder = radius .< 0.7
        @test all(radius[origin_cylinder] .≈ 0.6)
        @test sum(origin_cylinder) / nspins ≈ fraction_stuck / 3 rtol=0.1
        @test mean([mr.isinside(s.geometry, p) for p in stuck_spins[origin_cylinder]]) ≈ 0.5 rtol=0.1

        shifted_cylinder = radius .> 0.7
        other_radius = @. sqrt((abs(x) - 1)^2 + (abs(y) - 1)^2)
        @test all(other_radius[shifted_cylinder] .≈ 0.6)
        @test sum(shifted_cylinder) / nspins ≈ fraction_stuck * 2 / 3 rtol=0.1
        @test mean([mr.isinside(s.geometry, p) for p in stuck_spins[shifted_cylinder]]) ≈ 0.5 rtol=0.1
    end
end
@testset "Test simplified magnetisation transfer" begin
    @testset "Test signal loss between two walls" begin
        "Fraction of axons hitting the wall"
        frachit(x) = ((1 - erf(x)) + (1 - exp(-x*x)) / (sqrt(π) * x)) / 2
        frachit(wall_dist, diffusivity, timestep) = frachit(wall_dist / (2 * sqrt(diffusivity * timestep)))

        function test_MT_walls(wall_dist, diffusivity, timestep; nspins=100000, transfer=0.5)
            Random.seed!(1234)
            actual_transfer = 1 - (1 - transfer) ^ (1/sqrt(timestep))
            geometry = mr.walls(MT_fraction=actual_transfer)
            spins = [mr.Spin(position=Random.rand(3) .* wall_dist) for _ in 1:nspins]
            sequence = mr.Sequence(components=[mr.InstantRFPulse(flip_angle=90)], TR=1e5)
            simulation = mr.Simulation(sequence, geometry=geometry, diffusivity=diffusivity, max_timestep=timestep)
            signal = mr.signal(spins, simulation, [timestep])
            @test length(signal) == 1
            fhit = frachit(wall_dist, diffusivity, timestep)
            @test transfer * fhit ≈ (1 - mr.transverse(signal[end]) / nspins) rtol=0.1
        end
        test_MT_walls(2., 3., 0.1; transfer=0.5)
        test_MT_walls(2., 3., 0.1; transfer=0.1)
        test_MT_walls(10., 3., 0.1; transfer=0.5)
        test_MT_walls(10., 1., 0.1; transfer=0.5)
    end
    @testset "Test that transfer rate does not depend on timestep" begin
        Random.seed!(1234)
        geometry = mr.walls(repeats=1, MT_fraction=0.1)
        sequence = mr.Sequence(components=[mr.InstantRFPulse(flip_angle=90)], TR=1e5)

        reference = nothing
        for timestep in (0.01, 0.1, 1)
            simulation = mr.Simulation(sequence, geometry=geometry, diffusivity=1., max_timestep=timestep)
            signal = mean([mr.transverse(mr.evolve(10000, simulation, 10.)) for _ in 1:Int(timestep/0.01)]) / 10000
            if isnothing(reference)
                reference = signal
            else
                @test log(reference) ≈ log(signal) rtol=0.1
            end
        end
    end
end
@testset "Test realistic magnetisation transfer" begin
    @testset "Particles getting stuck reduces diffusivity" begin
        Random.seed!(123)
        geometry = mr.walls(repeats=1)
        for density in (0, 0.5, 1, 2)
            @testset "Density = $density" begin
                simulation = mr.Simulation([], geometry=geometry, diffusivity=1, surface_density=density, dwell_time=0.5)
                snap1 = mr.evolve(100000, simulation, 0)
                snap2 = mr.evolve(snap1, simulation, 10)
                displacement = mr.position.(snap2) .- mr.position.(snap1)

                fraction_stuck = density / (1 + density)
                @test var(map(d->d[2], displacement)) ≈ (1 - fraction_stuck) * 20 rtol=0.03
                @test var(map(d->d[3], displacement)) ≈ (1 - fraction_stuck) * 20 rtol=0.03
                @test all(map(d->(abs(d[1]) <= 1.00001), displacement))
                @test sum(mr.stuck.(snap1)) ≈ fraction_stuck * 1e5 rtol=0.03
                @test sum(mr.stuck.(snap2)) ≈ fraction_stuck * 1e5 rtol=0.03

                @test all(mr.isinside(simulation.geometry, snap1) .== mr.isinside(simulation.geometry, snap2))
            end
        end
    end
    @testset "Setting surface MRI properties" begin
        geometry = mr.walls(repeats=1, T2_surface=1e-6)
        init = mr.Snapshot(10000)
        seq = mr.Sequence(components=[mr.InstantRFPulse(flip_angle=90)], TR=1000)
        for density in (0, 0.5, 1)
            for dwell_time in (1, 2)
                @testset "Density = $density; dwell_time = $dwell_time" begin
                    simulation = mr.Simulation([seq], geometry=geometry, diffusivity=2, surface_density=density, dwell_time=dwell_time)
                    fraction_stuck = density / (1 + density)
                    release_rate = fraction_stuck / dwell_time
                    collision_rate = release_rate / (1 - fraction_stuck)
                    log_signal = log(mr.transverse(mr.evolve(init, simulation, 3)) / 1e4)
                    @test log_signal ≈ -collision_rate * 3 rtol=0.1
                end
            end
        end
    end
end
end