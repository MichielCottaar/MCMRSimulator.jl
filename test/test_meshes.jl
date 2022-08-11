@testset "Creating and using meshes" begin
    @testset "Computing normals" begin
        @test mr.normal(SA[0, 0, 0], SA[1, 0, 0], SA[0, 1, 0]) ≈ SA[0, 0, 1]
        @test mr.normal(SA[0, 0, 0], SA[2, 0, 0], SA[0, 1, 0]) ≈ SA[0, 0, 1]
        @test mr.normal(SA[0, 0, 0], SA[2, 1, 0], SA[-1, 1, 0]) ≈ SA[0, 0, 1]
        @test mr.normal(SA[0, 0, 1], SA[0, 1, 0], SA[1, 0, 0]) ≈ -SA[sqrt(1/3), sqrt(1/3), sqrt(1/3)]
        (p1, p2, p3) = SA[1, 4, 2], SA[0, 15.23, 3.1], SA[-213, 801.28380, 1]
        @test mr.normal(p1, p2, p3) ≈ mr.normal(p2, p3, p1)
        @test mr.normal(p1, p2, p3) ≈ mr.normal(p3, p1, p2)
        @test mr.normal(p1, p2, p3) ≈ -mr.normal(p1, p3, p2)
    end
    @testset "Bounding box calculation" begin
        bb = mr.BoundingBox(mr.box_mesh(grid_size=1))
        @test all(bb.lower < SA[-0.5, -0.5, -0.5])
        @test all(bb.upper > SA[0.5, 0.5, 0.5])
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
end