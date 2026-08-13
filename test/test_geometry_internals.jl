const GI = mr.Geometries.InternalNew
const PhysicalGeometry = GI.PhysicalGeometries.PhysicalGeometry
const GeometryVector = GI.PhysicalGeometries.Groups.GeometryVector
const GeometryTuple = GI.PhysicalGeometries.Groups.GeometryTuple
const Transformations = GI.PhysicalGeometries.Transformations
const Shift = Transformations.Shift
const Scale = Transformations.Scale
const Rotate = Transformations.Rotate
const Project = Transformations.Project
const BoundingBoxes = GI.InternalBoundingBoxes
const ObstructionIndex = GI.Indices.ObstructionIndex
const BaseObstructions = GI.PhysicalGeometries.BaseObstructions

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

    geometry_3d = TestGeometry{3}(0)
    geometry_2d = TestGeometry{2}(0)
    shift = Shift(geometry_3d, [1.0, 2.0, 3.0])
    position = SVector(4.0, 5.0, 6.0)
    normal = SVector(0.0, 1.0, 0.0)
    @test shift isa Transformations.Transformation{3, 3}
    @test Transformations.forward(shift, position) == SVector(5.0, 7.0, 9.0)
    @test Transformations.backward(shift, position) == SVector(3.0, 3.0, 3.0)
    @test Transformations.forward_normal(shift, normal) == normal
    @test Transformations.backward_normal(shift, normal) == normal

    intersection = GI.PhysicalGeometries.Intersection(
        0.5,
        normal,
        true,
        ObstructionIndex(SVector(1, 2)),
        false,
    )
    @test intersection.obstruction_index == ObstructionIndex(SVector(1, 2))
    @test !Base.isempty(intersection)
    @test Base.isempty(GI.PhysicalGeometries.Intersection{3}())

    infinite_wall = BaseObstructions.InfiniteWall()
    @test infinite_wall isa PhysicalGeometry{1}
    @test BoundingBoxes.InternalBoundingBox(infinite_wall) isa BoundingBoxes.InternalBoundingBox{1}
    wall_hit = GI.PhysicalGeometries.detect_intersection(
        infinite_wall,
        SVector(-1.0),
        SVector(1.0),
        ObstructionIndex(SVector(4)),
    )
    @test wall_hit.distance == 0.5
    @test wall_hit.normal == SVector(-1.0)
    @test wall_hit.obstruction_index == ObstructionIndex(SVector(4))

    cylinder = BaseObstructions.InfiniteCylinder(2.0)
    sphere = BaseObstructions.Sphere(2.0)
    @test BoundingBoxes.InternalBoundingBox(cylinder) isa BoundingBoxes.InternalBoundingBox{2}
    @test BoundingBoxes.InternalBoundingBox(sphere) isa BoundingBoxes.InternalBoundingBox{3}
    sphere_hit = GI.PhysicalGeometries.detect_intersection(
        sphere,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
    )
    @test sphere_hit.distance ≈ 1 / 6
    @test sphere_hit.normal == SVector(-1.0, 0.0, 0.0)

    triangle = BaseObstructions.FullTriangle(
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
    )
    triangle_box = BoundingBoxes.InternalBoundingBox(triangle)
    @test triangle_box isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.lower(triangle_box) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.upper(triangle_box) == SVector(1.0, 1.0, 0.0)
    triangle_hit = GI.PhysicalGeometries.detect_intersection(
        triangle,
        SVector(0.25, 0.25, 1.0),
        SVector(0.25, 0.25, -1.0),
    )
    @test triangle_hit.distance == 0.5
    @test triangle_hit.normal == SVector(0.0, 0.0, 1.0)

    scale = Scale(geometry_3d, 2.0)
    @test scale isa Transformations.Transformation{3, 3}
    @test scale.scale == 2.0
    @test Transformations.forward(scale, position) == SVector(8.0, 10.0, 12.0)
    @test Transformations.backward(scale, SVector(8.0, 10.0, 12.0)) == position
    @test Transformations.forward_normal(scale, normal) == normal
    @test Transformations.backward_normal(scale, normal) == normal
    @test_throws ArgumentError Scale(geometry_3d, 0.0)
    @test_throws ArgumentError Scale(geometry_3d, -1.0)
    rotation = Rotate(geometry_2d, SMatrix{2, 2, Float64}([0.0 -1.0; 1.0 0.0]))
    position_2d = SVector(1.0, 2.0)
    @test Transformations.forward(rotation, position_2d) == SVector(-2.0, 1.0)
    @test Transformations.backward(rotation, SVector(-2.0, 1.0)) == position_2d
    @test Transformations.forward_normal(rotation, position_2d) == SVector(-2.0, 1.0)

    projection = Project{3, 2}(geometry_2d)
    @test projection isa Transformations.Transformation{3, 2}
    @test Transformations.forward(projection, SVector(1.0, 2.0, 3.0)) == SVector(1.0, 2.0)
    @test_throws ArgumentError Transformations.backward(projection, SVector(1.0, 2.0))
    @test_throws ArgumentError Transformations.forward_normal(projection, SVector(1.0, 2.0, 3.0))
    @test_throws ArgumentError Transformations.backward_normal(projection, SVector(1.0, 2.0))

    geometry_1d = TestGeometry{1}(0)
    projection_z = Project{3, 1}(geometry_1d)
    @test Transformations.forward(projection_z, SVector(1.0, 2.0, 3.0)) == SVector(3.0)
    @test_throws ArgumentError Project{2, 1}(geometry_1d)

    cube = BoundingBoxes.InternalBoundingBox(2.0, [1.0, 2.0, 3.0])
    centered_cube = BoundingBoxes.InternalBoundingBox{3}(2.0)
    rect = BoundingBoxes.InternalBoundingBox([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
    centered_rect = BoundingBoxes.InternalBoundingBox{3}([1.0, 2.0, 3.0])

    for box in (cube, centered_cube, rect, centered_rect)
        @test box isa BoundingBoxes.InternalBoundingBox{3}
        @test BoundingBoxes.isinside(box, BoundingBoxes.lower(box))
        @test BoundingBoxes.isinside(box, BoundingBoxes.upper(box))
        @test BoundingBoxes.could_intersect(box, BoundingBoxes.lower(box), BoundingBoxes.upper(box))
        @test BoundingBoxes.does_intersect(box, BoundingBoxes.lower(box), BoundingBoxes.upper(box))
    end
    @test BoundingBoxes.center(centered_cube) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.center(cube) == SVector(1.0, 2.0, 3.0)
    @test BoundingBoxes.center(centered_rect) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.center(rect) == SVector(1.0, 2.0, 3.0)
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

    @test BoundingBoxes.InternalBoundingBox{3}(2.0) isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.InternalBoundingBox(2.0, [1.0, 2.0, 3.0]) isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.InternalBoundingBox{3}([1.0, 2.0, 3.0]) isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.InternalBoundingBox([1.0, 2.0, 3.0], [4.0, 5.0, 6.0]) isa BoundingBoxes.InternalBoundingBox{3}

    displacement = SVector(4.0, 5.0, 6.0)
    shifted_centered_cube = BoundingBoxes.shift(centered_cube, displacement)
    shifted_cube = BoundingBoxes.shift(cube, displacement)
    shifted_centered_rect = BoundingBoxes.shift(centered_rect, displacement)
    shifted_rect = BoundingBoxes.shift(rect, displacement)
    scaled_centered_cube = Transformations.forward(scale, centered_cube)
    scaled_rect = Transformations.forward(scale, rect)
    @test scaled_centered_cube isa BoundingBoxes.InternalBoundingBox{3}
    @test scaled_rect isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.lower(scaled_centered_cube) == SVector(-4.0, -4.0, -4.0)
    @test BoundingBoxes.upper(scaled_centered_cube) == SVector(4.0, 4.0, 4.0)
    @test BoundingBoxes.lower(scaled_rect) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.upper(scaled_rect) == SVector(4.0, 8.0, 12.0)
    @test Transformations.forward(shift, centered_cube) isa BoundingBoxes.InternalBoundingBox{3}
    transformed_centered_cube = Transformations.forward(shift, centered_cube)
    recovered_centered_cube = Transformations.backward(shift, transformed_centered_cube)
    @test BoundingBoxes.lower(recovered_centered_cube) == BoundingBoxes.lower(centered_cube)
    @test BoundingBoxes.upper(recovered_centered_cube) == BoundingBoxes.upper(centered_cube)
    @test shifted_centered_cube isa BoundingBoxes.InternalBoundingBox{3}
    @test shifted_cube isa BoundingBoxes.InternalBoundingBox{3}
    @test shifted_centered_rect isa BoundingBoxes.InternalBoundingBox{3}
    @test shifted_rect isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.lower(shifted_centered_cube) == SVector(2.0, 3.0, 4.0)
    @test BoundingBoxes.lower(shifted_cube) == SVector(3.0, 5.0, 7.0)
    @test BoundingBoxes.lower(shifted_centered_rect) == SVector(3.0, 3.0, 3.0)
    @test BoundingBoxes.lower(shifted_rect) == SVector(4.0, 5.0, 6.0)
    @test typeof(BoundingBoxes.shift(centered_cube, [4.0, 5.0, 6.0])) === typeof(shifted_centered_cube)
    @test Transformations.forward(projection, centered_cube) isa BoundingBoxes.InternalBoundingBox{2}
    @test Transformations.forward(projection_z, centered_cube) isa BoundingBoxes.InternalBoundingBox{1}
    @test_throws ArgumentError Transformations.backward(projection, centered_cube)

    rotation_3d = Rotate(geometry_3d, SMatrix{3, 3, Float64}(I))
    rotated_box = Transformations.forward(rotation_3d, centered_cube)
    @test rotated_box isa BoundingBoxes.InternalBoundingBox{3}
    transformed_bounds = (BoundingBoxes.lower(centered_cube), BoundingBoxes.upper(centered_cube))
    @test all(BoundingBoxes.isinside(rotated_box, Transformations.forward(rotation_3d, point)) for point in transformed_bounds)
    @test Transformations.backward(rotation_3d, rotated_box) isa BoundingBoxes.InternalBoundingBox{3}
end
