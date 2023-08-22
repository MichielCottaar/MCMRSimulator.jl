@testset "test_permeability.jl" begin
    for set_global in (false, true)
        @testset "Test single sphere for set_global=$set_global" begin
            Random.seed!(1234)
            snapshot = mr.Snapshot(zeros(10000, 3))

            function new_pos(permeability; timestep=0.1)
                if set_global
                    sphere = mr.spheres(radius=1.)
                    simulation = mr.Simulation([], geometry=sphere, diffusivity=1., max_timestep=timestep, permeability=permeability)
                else
                    sphere = mr.spheres(radius=1., permeability=permeability)
                    simulation = mr.Simulation([], geometry=sphere, diffusivity=1., max_timestep=timestep)
                end
                mr.position.(mr.evolve(snapshot, simulation, 10.))
            end
            pos = new_pos(0)
            @test maximum(norm.(pos)) < 1.

            pos = new_pos(1)
            @test maximum(norm.(pos)) > 1.
            for dim in 1:3
                @test std([a[dim] for a in pos]) ≈ sqrt(20.) rtol=0.1
            end

            for permeability in [0.1, 0.5, 0.8]
                pos = new_pos(permeability, timestep=0.01)
                @test maximum(norm.(pos)) > 1.
                for dim in 1:3
                    @test std([a[dim] for a in pos]) < sqrt(20.)
                end

                pos2 = new_pos(permeability, timestep=0.1)
                for dim in 1:3
                    @test std([a[dim] for a in pos]) ≈ std([a[dim] for a in pos2]) rtol=0.1
                end
            end
        end
        @testset "Test repeating wall for set_global=$set_global" begin
            Random.seed!(1234)
            snapshot = mr.Snapshot(zeros(10000, 3))

            function new_pos(permeability; timestep=0.1)
                if set_global
                    walls = mr.walls(position=0.5, repeats=1)
                    simulation = mr.Simulation([], geometry=walls, diffusivity=10., max_timestep=timestep, permeability=permeability)
                else
                    walls = mr.walls(position=0.5, permeability=permeability, repeats=1)
                    simulation = mr.Simulation([], geometry=walls, diffusivity=10., max_timestep=timestep)
                end
                [p[1] for p in mr.position.(mr.evolve(snapshot, simulation, 10.))]
            end
            pos = new_pos(0)
            @test maximum(abs.(pos)) < 0.5

            pos = new_pos(1)
            @test maximum(abs.(pos)) > 0.5
            @test std(pos) ≈ sqrt(200.) rtol=0.1

            for permeability in [0.1, 0.5, 0.8]
                pos = new_pos(permeability, timestep=0.01)
                @test maximum(abs.(pos)) > 0.5
                @test std(pos) < sqrt(200.)

                pos2 = new_pos(permeability, timestep=0.1)
                @test std(pos) ≈ std(pos2) rtol=0.1
            end
        end
    end
end
