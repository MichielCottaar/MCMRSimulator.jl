@testset "test_meshes.jl: Creating and using meshes" begin
    @testset "Computing normals" begin
        @test mr.normal([0, 0, 0], [1, 0, 0], [0, 1, 0]) ≈ [0, 0, 1]
        @test mr.normal([0, 0, 0], [2, 0, 0], [0, 1, 0]) ≈ [0, 0, 1]
        @test mr.normal([0, 0, 0], [2, 1, 0], [-1, 1, 0]) ≈ [0, 0, 1]
        @test mr.normal([0, 0, 1], [0, 1, 0], [1, 0, 0]) ≈ -[sqrt(1/3), sqrt(1/3), sqrt(1/3)]
        (p1, p2, p3) = [1., 4, 2], [0, 15.23, 3.1], [-213, 801.28380, 1]
        @test mr.normal(p1, p2, p3) ≈ mr.normal(p2, p3, p1)
        @test mr.normal(p1, p2, p3) ≈ mr.normal(p3, p1, p2)
        @test mr.normal(p1, p2, p3) ≈ -mr.normal(p1, p3, p2)
    end
    @testset "Bounding box calculation" begin
        bb = mr.BoundingBox(mr.box_mesh(grid_size=1))
        @test all(bb.lower < [-0.5, -0.5, -0.5])
        @test all(bb.upper > [0.5, 0.5, 0.5])
    end
    @testset "1x1x1 grid mesh intersection" begin
        mesh = mr.box_mesh(grid_size=1)
        @test size(mesh.grid) == (1, 1, 1)
        @test size(mesh.vertices) == (8, )
        @test size(mesh.triangles) == (12, )
        @test length(mesh.grid[1, 1, 1]) == 12
        @test mesh.grid[1, 1, 1] == collect(1:12)
    end
    @testset "Actual grid mesh intersection" begin
        mesh = mr.box_mesh(grid_size=5)
        @test size(mesh.grid) == (5, 5, 5)
        @test mesh.grid[2, 2, 2] == Int[]
        @test mesh.grid[3, 4, 2] == Int[]
        @test mesh.grid[3, 3, 1] == Int[1, 2]
        @test mesh.grid[2, 4, 1] == Int[1, 2]
        @test mesh.grid[1, 1, 1] == Int[1, 5, 6, 9, 10]
        @test mesh.grid[1, 3, 1] == Int[1, 6]
        @test mesh.grid[3, 1, 1] == Int[1, 9]
    end
    @testset "grid mesh box is closed" begin
        mesh = mr.box_mesh().grid
        for edge in (
            mesh[1, :, :],
            mesh[end, :, :],
            mesh[:, 1, :],
            mesh[:, end, :],
            mesh[:, :, 1],
            mesh[:, :, end],
        )
            @test all(map(l->length(l) > 0, edge))
        end
    end
    @testset "Simple bounces against mesh box" begin
        function compare(ms1 :: AbstractVector, ms2 :: AbstractVector)
            @test length(ms1) == length(ms2)
            for (m1, m2) in zip(ms1, ms2)
                @test m1 ≈ m2 atol=1e-4 rtol=1e-3
            end
        end
        @testset "Bounce on outside of box" begin
            mesh = mr.box_mesh()
            res = correct_collisions(
                mr.Movement([0.1, 0.1, 1], [0.1, 0.1, -1]),
                mesh
            )
            compare(res, [[0.1, 0.1, 1], [0.1, 0.1, 0.5], [0.1, 0.1, 2]])
        end
        @testset "Miss the box" begin
            mesh = mr.box_mesh()
            res = correct_collisions(
                mr.Movement([0, 0, 1.1], [0, 1.1, 0]),
                mesh
            )
            compare(res, [[0, 0, 1.1], [0, 1.1, 0]])
        end
        @testset "Straight bounce within the box" begin
            mesh = mr.box_mesh()
            res = correct_collisions(
                mr.Movement([0, 0, 0], [0, 0, 4]),
                mesh
            )
            compare(res, [
                [0, 0, 0], [0, 0, 0.5],
                [0, 0, -0.5], [0, 0, 0.5],
                [0, 0, -0.5], [0, 0, 0],
            ])
        end
    end
    @testset "Stay within a box" begin
        mesh = mr.box_mesh()
        snap = mr.Snapshot(rand(100, 3) .- 0.5)

        sequence = mr.dwi(bval=2., gradient_duration=0)

        simulation = mr.Simulation([sequence]; geometry=mesh, diffusivity=3.)
        final = mr.evolve(snap, simulation, 20)
        @test all(mr.position.(snap) != mr.position.(final))
        @test all(map(spin -> all(abs.(mr.position(spin) .<= 0.5)), final.spins))
    end
end