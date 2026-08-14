const FixBaseGeometry = mr.Geometries.Fix.FixBaseGeometry
const BaseObstructions = mr.Geometries.InternalNew.PhysicalGeometries.BaseObstructions

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
end
