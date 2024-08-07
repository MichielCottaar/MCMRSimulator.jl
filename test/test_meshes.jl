@testset "test_meshes.jl: Creating and using meshes" begin

function box_mesh(;center=[0, 0, 0], size=[1, 1, 1], kwargs...)
    center = SVector{3}(center)
    size = SVector{3}(size)
    vertices = [
        [0., 0, 0],
        [1., 0, 0],
        [0., 1, 0],
        [1., 1, 0],
        [0., 0, 1],
        [1., 0, 1],
        [0., 1, 1],
        [1., 1, 1],
    ]
    triangles = [
        # -z layer
        [1, 2, 3],
        [2, 3, 4],
        # +z layer
        [5, 6, 7],
        [6, 7, 8],
        # -x layer
        [1, 7, 5],
        [1, 3, 7],
        # +x layer
        [2, 8, 6],
        [2, 4, 8],
        # -y layer
        [1, 2, 6],
        [1, 5, 6],
        # +y layer
        [3, 4, 8],
        [3, 7, 8],
    ]
    c2 = center .- (size .* 0.5)
    return mr.Mesh(; vertices=map(v->((v .* size) .+ c2), vertices), triangles=triangles, kwargs...)
end

@testset "Computing normals" begin
    normal = mr.Geometries.Internal.Obstructions.Triangles.normal
    @test normal([0, 0, 0], [1, 0, 0], [0, 1, 0]) ≈ [0, 0, 1]
    @test normal([0, 0, 0], [2, 0, 0], [0, 1, 0]) ≈ [0, 0, 1]
    @test normal([0, 0, 0], [2, 1, 0], [-1, 1, 0]) ≈ [0, 0, 1]
    @test normal([0, 0, 1], [0, 1, 0], [1, 0, 0]) ≈ -[sqrt(1/3), sqrt(1/3), sqrt(1/3)]
    (p1, p2, p3) = [1., 4, 2], [0, 15.23, 3.1], [-213, 801.28380, 1]
    @test normal(p1, p2, p3) ≈ normal(p2, p3, p1)
    @test normal(p1, p2, p3) ≈ normal(p3, p1, p2)
    @test normal(p1, p2, p3) ≈ -normal(p1, p3, p2)
end

@testset "Computing inside for mesh" begin
    for center in ([0, 0, 0], [-0.11, 0, 0.], [-0.11, -0.11, -0.11])
        for grid_resolution in (0.1, 0.5, Inf)
            for repeats in (nothing, [1.5, 1.5, 1.5])
                mesh = box_mesh(center=[-0.11, -0.11, -0.11], grid_resolution=grid_resolution, repeats=repeats)
                @test mr.isinside(mesh, [0, 0, 0]) == 1
                @test mr.isinside(mesh, [-0.5, 0, 0]) == 1
                @test mr.isinside(mesh, [-0.5, 0, -0.5]) == 1
                @test mr.isinside(mesh, [-0.5, 0, 0.5]) == 0
                @test mr.isinside(mesh, [0.5, 0, 0.5]) == 0
                @test mr.isinside(mesh, [0.7, 0, -0.7]) == 0
            end
        end
    end
end
if false
    @testset "1x1x1 grid mesh intersection" begin
        mesh = box_mesh(grid_size=1)
        @test size(mesh.grid) == (1, 1, 1)
        @test size(mesh.vertices) == (8, )
        @test length(mesh) == 12
        @test length(mesh.grid[1, 1, 1]) == 12
        @test mesh.grid[1, 1, 1] == collect(1:12)
    end
    @testset "Actual grid mesh intersection" begin
        mesh = box_mesh(grid_size=5)
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
        mesh = box_mesh().grid
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
end
@testset "Simple bounces against mesh box" begin
    function compare(ms1 :: AbstractVector, ms2 :: AbstractVector)
        @test length(ms1) == length(ms2)
        for (m1, m2) in zip(ms1, ms2)
            @test m1 ≈ m2 atol=1e-4 rtol=1e-3
        end
    end
    @testset "Bounce on outside of box" begin
        mesh = box_mesh()
        res = correct_collisions(
            [0.1, 0.1, 1], [0.1, 0.1, -1],
            mesh
        )
        compare(res, [[0.1, 0.1, 1], [0.1, 0.1, 0.5], [0.1, 0.1, 2]])
    end
    @testset "Miss the box" begin
        mesh = box_mesh()
        res = correct_collisions(
            [0, 0, 1.1], [0, 1.1, 0],
            mesh
        )
        compare(res, [[0, 0, 1.1], [0, 1.1, 0]])
    end
    @testset "Straight bounce within the box" begin
        mesh = box_mesh()
        res = correct_collisions(
            [0, 0, 0], [0, 0, 4],
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
    mesh = box_mesh()
    snap = mr.Snapshot(rand(100, 3) .- 0.5)

    sequence = DWI(TE=80., bval=2., gradient=(type=:instant, ))

    simulation = mr.Simulation([sequence]; geometry=mesh, diffusivity=3.)
    final = mr.evolve(snap, simulation, 20)
    @test all(mr.position.(snap) != mr.position.(final))
    @test all(map(spin -> all(abs.(mr.position(spin) .<= 0.5)), final.spins))
end

@testset "Stay inside/outside a cylinder" begin
    Random.seed!(1)
    bendy_cylinder = mr.BendyCylinder(control_point=[0, 0, 0], radius=1., repeats=[2., 2., 2.], closed=[0, 0, 1])
    snap = mr.Snapshot(rand(100, 3) .* 2 .- 1.)

    real_inside(spin) = norm(spin.position[1:2]) .< 1.

    starts_inside = mr.isinside(bendy_cylinder, snap)
    @test all(real_inside.(snap) .== starts_inside)

    new_snap = mr.evolve(snap, mr.Simulation([], geometry=bendy_cylinder), 10.)

    @test all(real_inside.(snap) .== starts_inside)
    @test all(mr.isinside(bendy_cylinder, snap) .== starts_inside)
end

@testset "Compare mesh cylinder with perfect one" begin
    bendy_cylinder = mr.BendyCylinder(control_point=[0, 0, 0], radius=1., repeats=[2., 2., 2.], closed=[1, 0, 0])
    cylinder = mr.Cylinders(position=[0, 0], radius=1., repeats=[2., 2.], rotation=:x)
    fmesh = mr.fix(bendy_cylinder)
    fcylinder = mr.fix(cylinder)
    @testset "Compare inside" begin
        Random.seed!(1)
        for _ in 1:1000
            position = Random.rand(3) .* 1000.
            @test length(mr.isinside(fmesh, position)) == length(mr.isinside(fcylinder, position))
        end
    end
end

end