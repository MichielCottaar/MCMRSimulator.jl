@testset "test_permeability.jl" begin
    for set_global in (false, true)
        @testset "Test single sphere for set_global=$set_global" begin
            Random.seed!(1234)
            snapshot = mr.Snapshot(zeros(100000, 3))

            function new_pos(permeability; kwargs...)
                if set_global
                    sphere = mr.Spheres(radius=1.)
                    simulation = mr.Simulation([], geometry=sphere, diffusivity=1., permeability=permeability; timestep=(size_scale=1e5, kwargs...))
                else
                    sphere = mr.Spheres(radius=1., permeability=permeability)
                    simulation = mr.Simulation([], geometry=sphere, diffusivity=1.; timestep=(size_scale=1e5, kwargs...))
                end
                mr.position.(mr.evolve(snapshot, simulation, 10.))
            end
            init_pos = new_pos(0)
            @test maximum(norm.(init_pos)) < 1.
            snapshot = mr.Snapshot(init_pos)

            pos = new_pos(Inf)
            @test maximum(norm.(pos)) > 1.
            for dim in 1:3
                @test std([a[dim] for a in (pos .- init_pos)]) ≈ sqrt(20.) rtol=0.1
            end

            for permeability in [0.1, 1.]
                pos = new_pos(permeability)
                @test maximum(norm.(pos)) > 1.
                for dim in 1:3
                    @test std([a[dim] for a in (pos .- init_pos)]) < sqrt(20.)
                end

                pos2 = new_pos(permeability, permeability=0.05)
                for dim in 1:3
                    @test std([a[dim] for a in (pos .- init_pos)]) ≈ std([a[dim] for a in (pos2 .- init_pos)]) rtol=0.1
                end
            end
        end
        @testset "Test repeating wall for set_global=$set_global" begin
            Random.seed!(1234)
            snapshot = mr.Snapshot(zeros(10000, 3))

            function new_pos(permeability; timestep=0.1)
                if set_global
                    walls = mr.Walls(position=0.5, repeats=1)
                    simulation = mr.Simulation([], geometry=walls, diffusivity=10., timestep=timestep, permeability=permeability)
                else
                    walls = mr.Walls(position=0.5, permeability=permeability, repeats=1)
                    simulation = mr.Simulation([], geometry=walls, diffusivity=10., timestep=timestep)
                end
                [p[1] for p in mr.position.(mr.evolve(snapshot, simulation, 10.))]
            end
            pos = new_pos(0)
            @test maximum(abs.(pos)) < 0.5

            pos = new_pos(Inf)
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
