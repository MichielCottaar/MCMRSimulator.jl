const GI = mr.Geometries.InternalNew
const PhysicalGeometry = GI.PhysicalGeometries.PhysicalGeometry
const GeometryVector = GI.PhysicalGeometries.Groups.GeometryVector
const GeometryTuple = GI.PhysicalGeometries.Groups.GeometryTuple
const Transformations = GI.PhysicalGeometries.Transformations
const Shift = Transformations.Shift
const Rotate = Transformations.Rotate
const Project = Transformations.Project
const BoundingBoxes = GI.InternalBoundingBoxes

@testset "Geometry internals" begin
    struct TestGeometry{N} <: PhysicalGeometry{N}
        value::Int
    end

    geometries = [TestGeometry{3}(1), TestGeometry{3}(2)]
    vector = GeometryVector{3}(geometries)
    @test vector isa GI.PhysicalGeometries.Groups.GroupGeometry{3}
    @test size(vector) == (2,)
    @test axes(vector) == (Base.OneTo(2),)
    @test length(vector) == 2
    @test firstindex(vector) == 1
    @test lastindex(vector) == 2
    @test vector[1].value == 1
    @test [g.value for g in vector] == [1, 2]
    @test [g.value for g in vector[1:2]] == [1, 2]

    vector[2] = TestGeometry{3}(3)
    @test vector[2].value == 3

    tuple = GeometryTuple{3}((TestGeometry{3}(4), TestGeometry{3}(5)))
    @test tuple isa GI.PhysicalGeometries.Groups.GroupGeometry{3}
    @test length(tuple) == 2
    @test firstindex(tuple) == 1
    @test lastindex(tuple) == 2
    @test tuple[1].value == 4
    @test [g.value for g in tuple] == [4, 5]
    @test Tuple(tuple) == (TestGeometry{3}(4), TestGeometry{3}(5))

    @test length(GeometryVector{3}(TestGeometry{3}[])) == 0
    @test length(GeometryTuple{3}(())) == 0
    @test_throws MethodError GeometryVector{3}([TestGeometry{2}(1)])

    shift = Shift([1.0, 2.0, 3.0])
    position = SVector(4.0, 5.0, 6.0)
    normal = SVector(0.0, 1.0, 0.0)
    @test shift isa Transformations.Transformation{3, 3}
    @test Transformations.forward_position(shift, position) == SVector(5.0, 7.0, 9.0)
    @test Transformations.backward_position(shift, position) == SVector(3.0, 3.0, 3.0)
    @test Transformations.forward_normal(shift, normal) == normal
    @test Transformations.backward_normal(shift, normal) == normal

    rotation = Rotate(SMatrix{2, 2, Float64}([0.0 -1.0; 1.0 0.0]))
    position_2d = SVector(1.0, 2.0)
    @test Transformations.forward_position(rotation, position_2d) == SVector(-2.0, 1.0)
    @test Transformations.backward_position(rotation, SVector(-2.0, 1.0)) == position_2d
    @test Transformations.forward_normal(rotation, position_2d) == SVector(-2.0, 1.0)

    projection = Project{3, 2}()
    @test projection isa Transformations.Transformation{3, 2}
    @test Transformations.forward_position(projection, SVector(1.0, 2.0, 3.0)) == SVector(1.0, 2.0)
    @test_throws ArgumentError Transformations.backward_position(projection, SVector(1.0, 2.0))
    @test_throws ArgumentError Transformations.forward_normal(projection, SVector(1.0, 2.0, 3.0))
    @test_throws ArgumentError Transformations.backward_normal(projection, SVector(1.0, 2.0))

    projection_z = Project{3, 1}()
    @test Transformations.forward_position(projection_z, SVector(1.0, 2.0, 3.0)) == SVector(3.0)
    @test_throws ArgumentError Project{2, 1}()

    cube = BoundingBoxes.InternalBoundingCube([1.0, 2.0, 3.0], 2.0)
    centered_cube = BoundingBoxes.InternalCenteredBoundingCube{3}(2.0)
    rect = BoundingBoxes.InternalBoundingRect([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
    centered_rect = BoundingBoxes.InternalCenteredBoundingRect{3}([1.0, 2.0, 3.0])

    for box in (cube, centered_cube, rect, centered_rect)
        @test box isa BoundingBoxes.InternalBoundingBox{3}
        @test length(BoundingBoxes.corners(box)) == 8
        @test all(BoundingBoxes.isinside(box, corner) for corner in BoundingBoxes.corners(box))
        @test BoundingBoxes.could_intersect(box, BoundingBoxes.lower(box), BoundingBoxes.upper(box))
        @test BoundingBoxes.does_intersect(box, BoundingBoxes.lower(box), BoundingBoxes.upper(box))
    end
    @test BoundingBoxes.lower(cube) == SVector(-1.0, 0.0, 1.0)
    @test BoundingBoxes.upper(cube) == SVector(3.0, 4.0, 5.0)
    @test BoundingBoxes.lower(centered_cube) == SVector(-2.0, -2.0, -2.0)
    @test BoundingBoxes.upper(centered_cube) == SVector(2.0, 2.0, 2.0)
    @test BoundingBoxes.lower(rect) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.upper(rect) == SVector(2.0, 4.0, 6.0)
    @test BoundingBoxes.lower(centered_rect) == SVector(-1.0, -2.0, -3.0)
    @test BoundingBoxes.upper(centered_rect) == SVector(1.0, 2.0, 3.0)
    @test !BoundingBoxes.could_intersect(centered_cube, SVector(3.0, 3.0, 3.0), SVector(4.0, 4.0, 4.0))
    @test BoundingBoxes.could_intersect(
        centered_cube,
        SVector(-3.0, 3.0, 0.0),
        SVector(3.0, -100.0, 0.0),
    )
    @test !BoundingBoxes.does_intersect(
        centered_cube,
        SVector(-3.0, 3.0, 0.0),
        SVector(3.0, -100.0, 0.0),
    )

    @test BoundingBoxes.InternalBoundingBox(2.0) isa BoundingBoxes.InternalCenteredBoundingCube{3}
    @test BoundingBoxes.InternalBoundingBox(2.0, [1.0, 2.0, 3.0]) isa BoundingBoxes.InternalBoundingCube{3}
    @test BoundingBoxes.InternalBoundingBox([1.0, 2.0, 3.0]) isa BoundingBoxes.InternalCenteredBoundingRect{3}
    @test BoundingBoxes.InternalBoundingBox([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) isa BoundingBoxes.InternalBoundingRect{3}
end
