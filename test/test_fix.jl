const FixBaseGeometry = mr.Geometries.Fix.FixBaseGeometry
const FixTransformations = mr.Geometries.Fix.FixTransformations
const BaseObstructions = mr.Geometries.InternalNew.PhysicalGeometries.BaseObstructions
const Groups = mr.Geometries.InternalNew.PhysicalGeometries.Groups
const Repeats = mr.Geometries.InternalNew.PhysicalGeometries.Repeats
const Transformations = mr.Geometries.InternalNew.PhysicalGeometries.Transformations

@testset "Fix base geometry" begin
    @test FixBaseGeometry.fix_base_geometry(mr.Walls(position=0.0)) == [BaseObstructions.InfiniteWall()]
    @test FixBaseGeometry.fix_base_geometry(mr.Cylinders(radius=[1.0, 2.0])) == [
        BaseObstructions.InfiniteCylinder(1.0),
        BaseObstructions.InfiniteCylinder(2.0),
    ]
    @test FixBaseGeometry.fix_base_geometry(mr.Spheres(radius=[1.0, 2.0])) == [
        BaseObstructions.Sphere(1.0),
        BaseObstructions.Sphere(2.0),
    ]

    inner, outer = FixBaseGeometry.fix_base_geometry(mr.Annuli(inner=1.0, outer=2.0))
    @test inner == [BaseObstructions.InfiniteCylinder(1.0)]
    @test outer == [BaseObstructions.InfiniteCylinder(2.0)]

    overlapping = FixBaseGeometry.fix_base_geometry(mr.Spheres(
        radius=[1.0, 1.0],
        position=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]],
        overlapping=true,
    ))
    @test length(overlapping[1].overlaps_with) == 1
    @test length(overlapping[2].overlaps_with) == 1

    positioned = FixTransformations.fix_transformations(
        mr.Cylinders(radius=[1.0, 1.0], position=[[0.0, 0.0], [2.0, 0.0]]),
        FixBaseGeometry.fix_base_geometry(mr.Cylinders(radius=[1.0, 1.0], position=[[0.0, 0.0], [2.0, 0.0]])),
    )
    @test positioned isa Transformations.Rotate{3, 2}
    @test positioned.geometry isa Groups.GeometryVectorBoundingBox
    @test positioned.geometry.geometries[2].shift == SVector(-2.0, 0.0)

    global_position_group = mr.Spheres(radius=1.0, position=[1.0, 2.0, 3.0])
    globally_positioned = FixTransformations.fix_transformations(
        global_position_group,
        FixBaseGeometry.fix_base_geometry(global_position_group),
    )
    @test globally_positioned isa Transformations.Shift
    @test globally_positioned.shift == SVector(-1.0, -2.0, -3.0)

    grid_group = mr.Spheres(radius=1.0, grid_resolution=1.0)
    gridded = FixTransformations.fix_transformations(
        grid_group,
        FixBaseGeometry.fix_base_geometry(grid_group),
    )
    @test gridded isa Groups.GeometryVectorGrid

    repeated_group = mr.Spheres(radius=1.0, repeats=[4.0, 4.0, 4.0])
    repeated = FixTransformations.fix_transformations(
        repeated_group,
        FixBaseGeometry.fix_base_geometry(repeated_group),
    )
    @test repeated isa Repeats.Repeat

    rotated_group = mr.Cylinders(radius=1.0, rotation=:x)
    rotated = FixTransformations.fix_transformations(
        rotated_group,
        FixBaseGeometry.fix_base_geometry(rotated_group),
    )
    @test rotated isa Transformations.Rotate{3, 2}

    annulus_group = mr.Annuli(inner=1.0, outer=2.0, position=[0.0, 1.0])
    annulus_transformed = FixTransformations.fix_transformations(
        annulus_group,
        FixBaseGeometry.fix_base_geometry(annulus_group),
    )
    @test annulus_transformed isa Transformations.Rotate{3, 2}
    @test annulus_transformed.geometry isa Transformations.Shift{2}
    @test annulus_transformed.geometry.geometry isa Groups.GeometryTuple
    @test annulus_transformed.geometry.shift == SVector(0.0, -1.0)
end
