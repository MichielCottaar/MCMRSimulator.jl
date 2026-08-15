const GI = mr.Geometries.Internal
const PhysicalGeometry = GI.PhysicalGeometries.PhysicalGeometry
const GeometryVector = GI.PhysicalGeometries.Groups.GeometryVector
const GeometryVectorBoundingBox = GI.PhysicalGeometries.Groups.GeometryVectorBoundingBox
const GeometryVectorGrid = GI.PhysicalGeometries.Groups.GeometryVectorGrid
const GeometryTuple = GI.PhysicalGeometries.Groups.GeometryTuple
const Transformations = GI.PhysicalGeometries.Transformations
const Shift = Transformations.Shift
const Scale = Transformations.Scale
const Rotate = Transformations.Rotate
const Repeats = GI.PhysicalGeometries.Repeats
const Repeat = Repeats.Repeat
const BoundingBoxes = GI.InternalBoundingBoxes
const PublicBoundingBoxes = mr.Geometries.BoundingBoxes
const ObstructionIndex = GI.Indices.ObstructionIndex
const BaseObstructions = GI.PhysicalGeometries.BaseObstructions
const Meshes = GI.PhysicalGeometries.Meshes
const Mesh = Meshes.Mesh
const Transparents = GI.PhysicalGeometries.Transparents
const Transparent = GI.PhysicalGeometries.Transparent
const SizeScaleOverride = GI.SizeScales.SizeScaleOverride
const size_scale = GI.size_scale
const has_inside = GI.PhysicalGeometries.has_inside
const inside_indices = GI.PhysicalGeometries.inside_indices
const surface_sampling = GI.surface_sampling
const Properties = GI.Properties

@testset "public bounding boxes" begin
    box = PublicBoundingBoxes.BoundingBox([-1, -2, -3], [1, 2, 3])
    @test PublicBoundingBoxes.lower(box) == SVector(-1.0, -2.0, -3.0)
    @test PublicBoundingBoxes.upper(box) == SVector(1.0, 2.0, 3.0)
    @test PublicBoundingBoxes.BoundingBox(2.0) == PublicBoundingBoxes.BoundingBox(-2 .* ones(3), 2 .* ones(3))
    @test PublicBoundingBoxes.BoundingBox([1, 2, 3], 2.0) == PublicBoundingBoxes.BoundingBox([-1, 0, 1], [3, 4, 5])
    @test BoundingBoxes.InternalBoundingBox(box) isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.lower(BoundingBoxes.InternalBoundingBox(box)) == PublicBoundingBoxes.lower(box)
    @test BoundingBoxes.upper(BoundingBoxes.InternalBoundingBox(box)) == PublicBoundingBoxes.upper(box)
    @test_throws DimensionMismatch PublicBoundingBoxes.BoundingBox([0, 0], [1, 1])
    @test_throws ArgumentError PublicBoundingBoxes.BoundingBox(-1.0)
end
const get_value = Properties.get_value
const all_property_values = Properties.all_property_values

@testset "reflections" begin
    Reflections = mr.Reflections
    index = ObstructionIndex(SVector{2, Int}(1, 2))
    collision = GI.Intersection(0.5, SVector(1.0, 0.0, 0.0), false, index, false)

    free = Reflections.Reflection(1.0)
    @test !Reflections.has_intersection(free)
    @test isempty(Reflections.previous_hit(free).obstruction_index.indices)

    reflection = Reflections.Reflection(
        collision,
        SVector(-1.0, 0.0, 0.0),
        1.0,
        0.0,
        0.0,
    )
    @test Reflections.has_intersection(reflection)
    @test Reflections.has_hit(reflection) == index
    @test Reflections.previous_hit(reflection).obstruction_index == index
    @test reflection.direction == SVector(1.0, 0.0, 0.0)

    permeable = Reflections.Reflection(
        collision,
        SVector(-1.0, 0.0, 0.0),
        1.0,
        0.0,
        0.0,
        true,
    )
    @test permeable.inside
    @test permeable.direction == SVector(-1.0, 0.0, 0.0)
end

@testset "Geometry internals" begin
    @test collect(all_property_values(Properties.GeometryLeafProperties(1.0))) == [1.0]
    @test collect(all_property_values(Properties.GeometryVectorProperties([1.0, 2.0]))) == [1.0, 2.0]
    @test collect(all_property_values(
        Properties.GeometryVectorProperties([1.0, 2.0]),
        Properties.GeometryVectorProperties([3.0, 4.0]),
    )) == [(1.0, 3.0), (2.0, 4.0)]
    @test collect(all_property_values(
        Properties.GeometryTupleProperties((1.0, 2.0)),
        Properties.GeometryLeafProperties(3.0),
    )) == [(1.0, 3.0), (2.0, 3.0)]
    nested_properties = Properties.GeometryVectorProperties([
        Properties.GeometryTupleProperties((1.0, 2.0)),
        Properties.GeometryTupleProperties((3.0, 4.0)),
    ])
    @test collect(all_property_values(nested_properties)) == [1.0, 2.0, 3.0, 4.0]
    @test_throws ArgumentError collect(all_property_values(
        Properties.GeometryVectorProperties([1.0, 2.0]),
        Properties.GeometryTupleProperties((3.0, 4.0)),
    ))

    fixed = mr.fix([
        mr.Walls(position=0.0, permeability_surface=Inf, surface_density=0.0, dwell_time_surface=2.0),
        mr.Walls(position=1.0, permeability_surface=3.0, surface_density=2.0, dwell_time_surface=4.0),
    ]; surface_relaxation=5.0)
    @test GI.max_permeability_non_inf(fixed) == 3.0
    @test GI.max_surface_relaxation(fixed) == 5.0
    @test GI.min_dwell_time(fixed) == 4.0
    @test GI.max_timestep_sticking(fixed, 1.0, 0.1) ==
        0.1 * log(exp(-sqrt(π) * 2.0 / 4.0 / 2))^(-2)
    empty_fixed = mr.fix([])
    @test GI.max_permeability_non_inf(empty_fixed) == 0.0
    @test GI.max_surface_relaxation(empty_fixed) == 0.0
    @test GI.min_dwell_time(empty_fixed) == Inf
    @test GI.max_timestep_sticking(empty_fixed, 1.0, 0.1) == Inf

    density = Properties.GeometryLeafProperties(1.0)
    empty_positions, empty_normals = surface_sampling(
        BaseObstructions.Sphere(1.0),
        Properties.GeometryLeafProperties(0.0),
        BoundingBoxes.InternalBoundingBox{3}(2.0),
        10.0,
    )
    @test isempty(empty_positions)
    @test isempty(empty_normals)

    sphere_positions, sphere_normals = surface_sampling(
        BaseObstructions.Sphere(1.0),
        density,
        BoundingBoxes.InternalBoundingBox{3}(2.0),
        1.0,
    )
    @test length(sphere_positions) == length(sphere_normals)
    @test all(norm(position) ≈ 1.0 for position in sphere_positions)
    @test all(norm(normal) ≈ 1.0 for normal in sphere_normals)

    projected_cylinder = Transformations.Rotate(
        BaseObstructions.InfiniteCylinder(1.0),
        SMatrix{2, 3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0]),
    )
    projected_positions, projected_normals = surface_sampling(
        projected_cylinder,
        density,
        BoundingBoxes.InternalBoundingBox{3}(2.0),
        10.0,
    )
    @test length(projected_positions) == length(projected_normals)
    @test all(all(-2.0 .<= position .<= 2.0) for position in projected_positions)
    @test all(all(-2.0 .<= normal .<= 2.0) for normal in projected_normals)
    @test all(position[1]^2 + position[2]^2 ≈ 1.0 for position in projected_positions)
    @test all(norm(normal) ≈ 1.0 for normal in projected_normals)

    fixed_sample_geometry = mr.fix(mr.Spheres(radius=1.0); density=1.0)
    fixed_sample_positions, fixed_sample_normals = GI.surface_sampling(
        fixed_sample_geometry,
        BoundingBoxes.InternalBoundingBox{3}(2.0),
        100.0,
    )
    @test length(fixed_sample_positions) == length(fixed_sample_normals)
    @test all(norm(position) ≈ 1.0 for position in fixed_sample_positions)
    @test all(norm(normal) ≈ 1.0 for normal in fixed_sample_normals)
    empty_sample_positions, empty_sample_normals = GI.surface_sampling(
        mr.fix([]),
        BoundingBoxes.InternalBoundingBox{3}(2.0),
        100.0,
    )
    @test isempty(empty_sample_positions)
    @test isempty(empty_sample_normals)

    property_geometry = mr.fix(
        mr.Spheres(
            radius=[1.0, 2.0],
            permeability_surface=[3.0, 4.0],
            density_surface=[5.0, 6.0],
            dwell_time_surface=[7.0, 8.0],
            relaxation_surface=[9.0, 10.0],
        ),
    )
    property_intersection = GI.PhysicalGeometries.Intersection(
        0.5,
        SVector(1.0, 0.0, 0.0),
        false,
        ObstructionIndex(SVector(2, 1)),
        false,
    )
    @test GI.permeability(property_geometry, property_intersection) == 4.0
    @test GI.surface_relaxation(property_geometry, property_intersection) == 10.0
    @test GI.surface_density(property_geometry, property_intersection) == 6.0
    @test GI.dwell_time(property_geometry, property_intersection) == 8.0
    gap_intersection = GI.PhysicalGeometries.Intersection(
        property_intersection.distance,
        property_intersection.normal,
        property_intersection.inside,
        property_intersection.obstruction_index,
        true,
    )
    @test GI.permeability(property_geometry, gap_intersection) == Inf
    @test GI.surface_relaxation(property_geometry, gap_intersection) == 0.0
    @test GI.surface_density(property_geometry, gap_intersection) == 0.0
    @test GI.dwell_time(property_geometry, gap_intersection) == 0.0
    empty_intersection = GI.PhysicalGeometries.Intersection{3}()
    @test_throws AssertionError GI.permeability(property_geometry, empty_intersection)
    @test_throws AssertionError GI.surface_relaxation(property_geometry, empty_intersection)
    @test_throws AssertionError GI.surface_density(property_geometry, empty_intersection)
    @test_throws AssertionError GI.dwell_time(property_geometry, empty_intersection)

    mri_geometry = mr.fix(mr.Spheres(
        radius=[1.0, 2.0],
        R1_inside=[1.0, 2.0],
        R2_inside=[3.0, 4.0],
        off_resonance_inside=[5.0, 6.0],
        R1_surface=[10.0, 20.0],
        R2_surface=[30.0, 40.0],
        off_resonance_surface=[50.0, 60.0],
    ))
    inside_second = GI.IsInside([ObstructionIndex(SVector(2))])
    mri_intersection = GI.PhysicalGeometries.Intersection(
        0.5,
        SVector(1.0, 0.0, 0.0),
        false,
        ObstructionIndex(SVector(2)),
        false,
    )
    @test GI.R1(mri_geometry, inside_second) == 2.0
    @test GI.R2(mri_geometry, inside_second) == 4.0
    @test GI.off_resonance(mri_geometry, inside_second) == 6.0
    @test GI.R1(mri_geometry, inside_second, mri_intersection) == 22.0
    @test GI.R2(mri_geometry, inside_second, mri_intersection) == 44.0
    @test GI.off_resonance(mri_geometry, inside_second, mri_intersection) == 66.0

    scalar_mri_geometry = mr.fix(mr.Spheres(
        radius=[1.0, 2.0],
        R1_inside=3.0,
        R2_inside=4.0,
        off_resonance_inside=5.0,
    ))
    inside_both = GI.IsInside([
        ObstructionIndex(SVector(1)),
        ObstructionIndex(SVector(2)),
    ])
    @test GI.R1(scalar_mri_geometry, inside_both) == 3.0
    @test GI.R2(scalar_mri_geometry, inside_both) == 4.0
    @test GI.off_resonance(scalar_mri_geometry, inside_both) == 5.0

    @test GI.R1(mri_geometry, SVector(0.0, 0.0, 0.0)) == 3.0

    struct TestGeometry{N} <: PhysicalGeometry{N}
        value::Int
    end

    geometries = [TestGeometry{3}(1), TestGeometry{3}(2)]
    vector = GeometryVector{3}(geometries)
    @test GeometryVector(geometries) isa GeometryVector{3, TestGeometry{3}}
    @test GeometryVector([BaseObstructions.Sphere(1.0)]; bounding_box=true) isa GeometryVectorBoundingBox{3, BaseObstructions.Sphere}
    @test GeometryVector([BaseObstructions.Sphere(1.0)]; grid=true) isa GeometryVectorGrid{3, BaseObstructions.Sphere}
    @test_throws ArgumentError GeometryVector([BaseObstructions.Sphere(1.0)]; bounding_box=true, grid=true)
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
    @test has_inside(typeof(BaseObstructions.Sphere(1.0)))
    @test has_inside(typeof(BaseObstructions.OverlappingSphere(1.0)))
    @test !has_inside(typeof(BaseObstructions.InfiniteWall()))
    @test !has_inside(typeof(BaseObstructions.FullTriangle(
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
    )))
    mesh = Mesh(
        [
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0),
            SVector(0.0, 0.0, 1.0),
        ],
        [
            SVector(1, 2, 3),
            SVector(1, 2, 4),
            SVector(1, 3, 4),
        ],
        grid_resolution=0.5,
    )
    @test mesh isa PhysicalGeometry{3}
    @test has_inside(typeof(mesh))
    @test mesh.vertices[2] == SVector(1.0, 0.0, 0.0)
    @test mesh.first_index_of_gap == 4
    @test length(mesh.indices) == 4
    render_mesh = GI.geometry_mesh(mesh)
    @test length(render_mesh) == 1
    @test render_mesh[1].vertices == mesh.vertices
    @test render_mesh[1].triangles == mesh.indices[1:(mesh.first_index_of_gap - 1)]
    shifted_render_mesh = GI.geometry_mesh(Shift(mesh, [1.0, 2.0, 3.0]))
    @test shifted_render_mesh[1].vertices == [vertex - SVector(1.0, 2.0, 3.0) for vertex in mesh.vertices]
    grouped_render_mesh = GI.geometry_mesh(GeometryTuple{3}((mesh, mesh)))
    @test length(grouped_render_mesh) == 2
    repeated_mesh = Repeat(mesh, [2.0, 2.0, 2.0])
    @test length(GI.geometry_mesh(repeated_mesh)) == 1
    bounded_repeated_mesh = GI.geometry_mesh(
        repeated_mesh,
        PublicBoundingBoxes.BoundingBox([1.0, 0.0, 0.0], [2.0, 1.0, 1.0]),
    )
    @test length(bounded_repeated_mesh) == 2
    @test sort([minimum(vertex[1] for vertex in render_mesh.vertices) for render_mesh in bounded_repeated_mesh]) == [0.0, 2.0]
    @test all(length(render_mesh.triangles) == mesh.first_index_of_gap - 1 for render_mesh in bounded_repeated_mesh)
    transformed_bounded_mesh = GI.geometry_mesh(
        Shift(repeated_mesh, [1.0, 0.0, 0.0]),
        PublicBoundingBoxes.BoundingBox([-1.0, 0.0, 0.0], [0.0, 1.0, 1.0]),
    )
    @test length(transformed_bounded_mesh) == 1
    @test minimum(vertex[1] for vertex in transformed_bounded_mesh[1].vertices) == -1.0
    @test inside_indices(mesh, SVector(0.2, 0.2, 0.2)) == [ObstructionIndex()]
    @test inside_indices(mesh, SVector(1.2, 0.2, 0.2)) == ObstructionIndex[]
    @test mesh.indices[mesh.first_index_of_gap] == SVector(4, 2, 3)
    @test BoundingBoxes.lower(mesh.bounding_box) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.upper(mesh.bounding_box) == SVector(1.0, 1.0, 1.0)
    @test size_scale(mesh) == 1 / (2 * Meshes.curvature(mesh))
    @test size_scale(SizeScaleOverride(mesh, 0.25)) == 0.25
    @test SizeScaleOverride(mesh, 0.25) isa Transparent
    @test size_scale(Transformations.Scale(mesh, 2.0)) == 2 * size_scale(mesh)
    @test size_scale(Repeats.Repeat(mesh, [4.0, 4.0, 4.0])) == size_scale(mesh)
    @test Meshes.triangle(mesh, 1) == BaseObstructions.FullTriangle(
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
    )
    @test Meshes.triangle(mesh, 3) == BaseObstructions.FullTriangle(
        SVector(0.0, 1.0, 0.0),
        SVector(0.0, 0.0, 0.0),
        SVector(0.0, 0.0, 1.0),
    )
    @test Meshes.curvature(mesh) == Meshes.curvature(
        mesh.indices[1:(mesh.first_index_of_gap - 1)],
        mesh.vertices,
    )
    @test Meshes.curvature(mesh; include_gap_triangles=true) == Meshes.curvature(
        mesh.indices,
        mesh.vertices,
    )
    mesh_real_hit = GI.PhysicalGeometries.detect_intersection(
        mesh,
        SVector(0.2, 0.2, -1.0),
        SVector(0.2, 0.2, 0.2),
    )
    @test mesh_real_hit.distance ≈ 5 / 6
    @test mesh_real_hit.obstruction_index == ObstructionIndex(SVector(1))
    @test !mesh_real_hit.hit_gap
    mesh_gap_hit = GI.PhysicalGeometries.detect_intersection(
        mesh,
        SVector(0.2, 0.2, 0.2),
        SVector(1.0, 1.0, 1.0),
    )
    @test mesh_gap_hit.distance ≈ 1 / 6
    @test mesh_gap_hit.obstruction_index == ObstructionIndex(SVector(4))
    @test mesh_gap_hit.hit_gap
    @test Base.isempty(GI.PhysicalGeometries.detect_intersection(
        mesh,
        SVector(0.2, 0.2, 0.2),
        SVector(1.0, 1.0, 1.0),
        mesh_gap_hit,
    ))
    @test_throws ArgumentError Mesh(
        [SVector(0.0, 0.0, 0.0)],
        [SVector(1, 2, 3)],
    )
    @test has_inside(typeof(Shift(BaseObstructions.Sphere(1.0), [0.0, 0.0, 0.0])))
    @test !has_inside(typeof(Scale(BaseObstructions.InfiniteWall(), 1.0)))
    @test has_inside(typeof(GeometryVector{3}([BaseObstructions.Sphere(1.0)])))
    @test !has_inside(typeof(GeometryVector{1}([BaseObstructions.InfiniteWall()])))
    @test has_inside(typeof(GeometryTuple{3}((
        BaseObstructions.FullTriangle(
            SVector(0.0, 0.0, 0.0),
            SVector(1.0, 0.0, 0.0),
            SVector(0.0, 1.0, 0.0),
        ),
        BaseObstructions.Sphere(1.0),
    ))))
    @test !has_inside(typeof(GeometryTuple{3}((BaseObstructions.FullTriangle(
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
    ),))))
    @test !has_inside(typeof(GeometryTuple{3}(())))
    repeated_sphere = Repeat(BaseObstructions.Sphere(1.0), [4.0, 4.0, 4.0])
    @test has_inside(typeof(repeated_sphere))
    @test repeated_sphere.lower_overlap == SVector(0.0, 0.0, 0.0)
    @test repeated_sphere.upper_overlap == SVector(0.0, 0.0, 0.0)
    @test inside_indices(repeated_sphere, SVector(4.5, 0.0, 0.0)) == [ObstructionIndex()]
    @test inside_indices(repeated_sphere, SVector(2.0, 0.0, 0.0)) == ObstructionIndex[]
    @test_throws ArgumentError Repeat(BaseObstructions.Sphere(1.0), [0.0, 4.0, 4.0])
    @test_throws ArgumentError Repeat(BaseObstructions.Sphere(2.0), [1.0, 4.0, 4.0])
    @test_throws ArgumentError Repeat(
        Shift(BaseObstructions.Sphere(1.0), [3.0, 0.0, 0.0]),
        [2.0, 4.0, 4.0],
    )
    @test_throws ArgumentError BoundingBoxes.InternalBoundingBox(repeated_sphere)
    repeated_hit = GI.PhysicalGeometries.detect_intersection(
        repeated_sphere,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    @test repeated_hit.distance ≈ 0.4
    @test repeated_hit.normal ≈ SVector(-1.0, 0.0, 0.0)
    overlapping_repeat = Repeat(BaseObstructions.Sphere(1.5), [2.0, 4.0, 4.0])
    @test overlapping_repeat.lower_overlap == SVector(0.5, 0.0, 0.0)
    @test overlapping_repeat.upper_overlap == SVector(0.5, 0.0, 0.0)
    @test inside_indices(overlapping_repeat, SVector(1.4, 0.0, 0.0)) == [ObstructionIndex()]
    overlapping_hit = GI.PhysicalGeometries.detect_intersection(
        overlapping_repeat,
        SVector(1.4, 0.0, 0.0),
        SVector(1.9, 0.0, 0.0),
    )
    @test overlapping_hit.distance ≈ 0.2

    @test get_value(3.0, ObstructionIndex(SVector(1, 2))) == 3.0
    leaf_properties = Properties.GeometryLeafProperties(3.0)
    @test leaf_properties isa Properties.GeometryProperties{Float64}
    @test eltype(leaf_properties) === Float64
    vector_properties = Properties.GeometryVectorProperties([10.0, 20.0])
    @test vector_properties isa Properties.GeometryProperties{Float64}
    @test eltype(vector_properties) === Float64
    @test get_value(vector_properties, ObstructionIndex(SVector(2, 1))) == 20.0
    tuple_properties = Properties.GeometryTupleProperties((1.0, vector_properties))
    @test tuple_properties isa Properties.GeometryProperties{Float64}
    @test eltype(tuple_properties) === Float64
    @test get_value(tuple_properties, ObstructionIndex(SVector(2, 2))) == 20.0
    @test get_value(tuple_properties, [
        ObstructionIndex(SVector(1, 4)),
        ObstructionIndex(SVector(2, 2)),
    ]) == 21.0
    aggregation_properties = Properties.GeometryTupleProperties((
        Properties.GeometryVectorProperties([1.0, 2.0]),
        10.0,
    ))
    @test get_value(aggregation_properties, [
        ObstructionIndex(SVector(1, 1)),
        ObstructionIndex(SVector(1, 2)),
        ObstructionIndex(SVector(2, 1)),
        ObstructionIndex(SVector(2, 2)),
    ]) == 13.0
    @test get_value(aggregation_properties, ObstructionIndex[]) == 0
    @test_throws ArgumentError get_value(aggregation_properties, [
        ObstructionIndex(SVector(2, 1)),
        ObstructionIndex(SVector(1, 1)),
    ])
    @test_throws BoundsError get_value(vector_properties, ObstructionIndex(SVector(3)))
    @test_throws ArgumentError Properties.GeometryVectorProperties([
        Properties.GeometryLeafProperties(1.0),
        Properties.GeometryLeafProperties(2),
    ])
    @test_throws ArgumentError Properties.GeometryTupleProperties((
        Properties.GeometryLeafProperties(1.0),
        Properties.GeometryLeafProperties(2),
    ))

    sphere = BaseObstructions.Sphere(1.0)
    @test inside_indices(sphere, SVector(0.0, 0.0, 0.0)) == [ObstructionIndex()]
    @test inside_indices(sphere, SVector(2.0, 0.0, 0.0)) == ObstructionIndex[]
    @test inside_indices(BaseObstructions.InfiniteWall(), SVector(0.0)) == ObstructionIndex[]

    shifted_inside = Shift(sphere, [2.0, 0.0, 0.0])
    @test inside_indices(shifted_inside, SVector(-2.0, 0.0, 0.0)) == [ObstructionIndex()]
    @test inside_indices(shifted_inside, SVector(0.0, 0.0, 0.0)) == ObstructionIndex[]

    inside_spheres = GeometryVector{3}([
        Shift(sphere, [-0.5, 0.0, 0.0]),
        Shift(sphere, [0.5, 0.0, 0.0]),
    ])
    @test inside_indices(inside_spheres, SVector(0.0, 0.0, 0.0)) == [
        ObstructionIndex(SVector(1)),
        ObstructionIndex(SVector(2)),
    ]
    walls = GeometryVector{1}([BaseObstructions.InfiniteWall(), BaseObstructions.InfiniteWall()])
    @test inside_indices(walls, SVector(0.0)) == ObstructionIndex[]

    nested = GeometryTuple{3}((inside_spheres, sphere))
    @test inside_indices(nested, SVector(0.0, 0.0, 0.0)) == [
        ObstructionIndex(SVector(1, 1)),
        ObstructionIndex(SVector(1, 2)),
        ObstructionIndex(SVector(2)),
    ]
    @test_throws MethodError GeometryVector{3}([TestGeometry{2}(1)])
    @test_throws ArgumentError GeometryVector{3, PhysicalGeometry{3}}(PhysicalGeometry{3}[])
    @test_throws ArgumentError GeometryVector{3, Shift{3, PhysicalGeometry{3}}}(
        Shift{3, PhysicalGeometry{3}}[],
    )

    geometry_3d = TestGeometry{3}(0)
    geometry_2d = TestGeometry{2}(0)
    shift = Shift(geometry_3d, [1.0, 2.0, 3.0])
    position = SVector(4.0, 5.0, 6.0)
    normal = SVector(0.0, 1.0, 0.0)
    @test shift isa Transformations.Transformation{3, 3, typeof(geometry_3d)}
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
    )
    @test wall_hit.distance == 0.5
    @test wall_hit.normal == SVector(-1.0)
    previous_wall_hit = GI.PhysicalGeometries.Intersection(
        0.5,
        SVector(1.0, 0.0, 0.0),
        false,
        ObstructionIndex(),
        false,
    )
    @test Base.isempty(GI.PhysicalGeometries.detect_intersection(
        infinite_wall,
        SVector(-1.0),
        SVector(1.0),
        previous_wall_hit,
    ))

    cylinder = BaseObstructions.InfiniteCylinder(2.0)
    sphere = BaseObstructions.Sphere(2.0)
    @test BoundingBoxes.InternalBoundingBox(cylinder) isa BoundingBoxes.InternalBoundingBox{2}
    @test BoundingBoxes.InternalBoundingBox(sphere) isa BoundingBoxes.InternalBoundingBox{3}
    sphere_box = BoundingBoxes.InternalBoundingBox(sphere)
    @test BoundingBoxes.lower(sphere_box) == SVector(-2.0, -2.0, -2.0)
    @test BoundingBoxes.upper(sphere_box) == SVector(2.0, 2.0, 2.0)
    sphere_hit = GI.PhysicalGeometries.detect_intersection(
        sphere,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
    )
    @test sphere_hit.distance ≈ 1 / 6
    @test sphere_hit.normal == SVector(-1.0, 0.0, 0.0)
    boundary = SVector(2.0, 0.0, 0.0)
    boundary_inside = GI.PhysicalGeometries.Intersection(
        0.5,
        SVector(1.0, 0.0, 0.0),
        true,
        ObstructionIndex(),
        false,
    )
    boundary_outside = GI.PhysicalGeometries.Intersection(
        0.5,
        SVector(1.0, 0.0, 0.0),
        false,
        ObstructionIndex(),
        false,
    )
    @test inside_indices(sphere, boundary, boundary_inside) == [ObstructionIndex()]
    @test inside_indices(sphere, boundary, boundary_outside) == ObstructionIndex[]

    fixed_sphere = mr.fix(mr.Spheres(radius=2.0))
    fixed_boundary_inside = GI.PhysicalGeometries.Intersection(
        boundary_inside.distance,
        boundary_inside.normal,
        boundary_inside.inside,
        ObstructionIndex(SVector(1)),
        boundary_inside.hit_gap,
    )
    fixed_boundary_outside = GI.PhysicalGeometries.Intersection(
        boundary_outside.distance,
        boundary_outside.normal,
        boundary_outside.inside,
        ObstructionIndex(SVector(1)),
        boundary_outside.hit_gap,
    )
    @test GI.isinside(fixed_sphere, boundary, fixed_boundary_inside).inside_of == [ObstructionIndex(SVector(1))]
    @test GI.isinside(fixed_sphere, boundary, fixed_boundary_outside).inside_of == ObstructionIndex[]

    @test Base.isempty(GI.PhysicalGeometries.detect_intersection(
        sphere,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
        sphere_hit,
    ))
    previous_inside_sphere = GI.PhysicalGeometries.Intersection(
        0.5,
        SVector(1.0, 0.0, 0.0),
        true,
        ObstructionIndex(),
        false,
    )
    inside_sphere_hit = GI.PhysicalGeometries.detect_intersection(
        sphere,
        SVector(-1.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
        previous_inside_sphere,
    )
    @test inside_sphere_hit.distance ≈ 3 / 4
    @test inside_sphere_hit.normal == SVector(-1.0, 0.0, 0.0)

    fixed_sphere_hit = GI.detect_intersection(
        fixed_sphere,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
    )
    underlying_sphere_hit = GI.PhysicalGeometries.detect_intersection(
        fixed_sphere.geometry,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
    )
    @test fixed_sphere_hit == underlying_sphere_hit

    shifted_sphere = Shift(sphere, [1.0, 0.0, 0.0])
    shifted_sphere_box = BoundingBoxes.InternalBoundingBox(shifted_sphere)
    @test BoundingBoxes.lower(shifted_sphere_box) == SVector(-3.0, -2.0, -2.0)
    @test BoundingBoxes.upper(shifted_sphere_box) == SVector(1.0, 2.0, 2.0)
    shifted_hit = GI.PhysicalGeometries.detect_intersection(
        shifted_sphere,
        SVector(-4.0, 0.0, 0.0),
        SVector(2.0, 0.0, 0.0),
    )
    @test shifted_hit.distance ≈ 1 / 6
    @test shifted_hit.normal == SVector(-1.0, 0.0, 0.0)

    grouped_spheres = GeometryVector{3}([
        Shift(BaseObstructions.Sphere(1.0), [-2.0, 0.0, 0.0]),
        Shift(BaseObstructions.Sphere(1.0), [2.0, 0.0, 0.0]),
    ])
    grouped_box = BoundingBoxes.InternalBoundingBox(grouped_spheres)
    @test BoundingBoxes.lower(grouped_box) == SVector(-3.0, -1.0, -1.0)
    @test BoundingBoxes.upper(grouped_box) == SVector(3.0, 1.0, 1.0)
    filtered_spheres = GeometryVectorBoundingBox([
        Shift(sphere, [-0.5, 0.0, 0.0]),
        Shift(sphere, [5.0, 0.0, 0.0]),
    ])
    @test filtered_spheres isa GI.PhysicalGeometries.Groups.GeometryVectorLike{3}
    @test length(filtered_spheres.bounding_boxes) == 2
    filtered_hit = GI.PhysicalGeometries.detect_intersection(
        filtered_spheres,
        SVector(-2.0, 0.0, 0.0),
        SVector(2.0, 0.0, 0.0),
    )
    @test filtered_hit.obstruction_index == ObstructionIndex(SVector(1))
    @test inside_indices(filtered_spheres, SVector(0.0, 0.0, 0.0)) == [
        ObstructionIndex(SVector(1)),
    ]
    grid_spheres = GeometryVector([
        Shift(sphere, [-3.0, 0.0, 0.0]),
        Shift(sphere, [3.0, 0.0, 0.0]),
    ]; grid=true, grid_resolution=2.0)
    @test grid_spheres isa GeometryVectorGrid{3, Shift{3, BaseObstructions.Sphere}}
    @test inside_indices(grid_spheres, SVector(3.0, 0.0, 0.0)) == [ObstructionIndex(SVector(1))]
    regular_grid_hit = GI.PhysicalGeometries.detect_intersection(
        GeometryVector([
            Shift(sphere, [-3.0, 0.0, 0.0]),
            Shift(sphere, [3.0, 0.0, 0.0]),
        ]),
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    grid_hit = GI.PhysicalGeometries.detect_intersection(
        grid_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    @test grid_hit.distance == regular_grid_hit.distance
    @test grid_hit.obstruction_index == regular_grid_hit.obstruction_index
    @test_throws ArgumentError GeometryVectorBoundingBox(
        [sphere],
        BoundingBoxes.InternalBoundingBox{3}[],
    )
    @test_throws ArgumentError BoundingBoxes.InternalBoundingBox(GeometryVector{3}(TestGeometry{3}[]))
    @test_throws ArgumentError BoundingBoxes.InternalBoundingBox(GeometryTuple{3}(()))
    grouped_hit = GI.PhysicalGeometries.detect_intersection(
        grouped_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    @test grouped_hit.distance ≈ 1 / 5
    @test grouped_hit.obstruction_index == ObstructionIndex(SVector(2))
    grouped_boundary = SVector(3.0, 0.0, 0.0)
    grouped_inside = GI.PhysicalGeometries.Intersection(
        0.5,
        SVector(1.0, 0.0, 0.0),
        true,
        ObstructionIndex(SVector(2)),
        false,
    )
    grouped_outside = GI.PhysicalGeometries.Intersection(
        0.5,
        SVector(1.0, 0.0, 0.0),
        false,
        ObstructionIndex(SVector(2)),
        false,
    )
    @test inside_indices(grouped_spheres, grouped_boundary, grouped_inside) == [ObstructionIndex(SVector(2))]
    @test inside_indices(grouped_spheres, grouped_boundary, grouped_outside) == ObstructionIndex[]

    next_grouped_hit = GI.PhysicalGeometries.detect_intersection(
        grouped_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
        grouped_hit,
    )
    @test next_grouped_hit.distance ≈ 3 / 5
    @test next_grouped_hit.obstruction_index == ObstructionIndex(SVector(1))
    @test_throws ArgumentError GI.PhysicalGeometries.detect_intersection(
        grouped_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
        GI.PhysicalGeometries.Intersection(0.5, normal, false, ObstructionIndex(), false),
    )
    @test_throws ArgumentError GI.PhysicalGeometries.detect_intersection(
        grouped_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
        GI.PhysicalGeometries.Intersection(0.5, normal, false, ObstructionIndex(SVector(3)), false),
    )

    scaled_sphere = Scale(sphere, 2.0)
    scaled_hit = GI.PhysicalGeometries.detect_intersection(
        scaled_sphere,
        SVector(-6.0, 0.0, 0.0),
        SVector(6.0, 0.0, 0.0),
    )
    @test scaled_hit.distance ≈ 5 / 12
    @test scaled_hit.normal ≈ SVector(-1.0, 0.0, 0.0)

    projected_cylinder = Rotate(cylinder, SMatrix{2, 3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0]))
    projected_hit = GI.PhysicalGeometries.detect_intersection(
        projected_cylinder,
        SVector(-3.0, 0.0, 4.0),
        SVector(3.0, 0.0, 4.0),
    )
    @test projected_hit.distance ≈ 1 / 6
    @test projected_hit.normal == SVector(-1.0, 0.0, 0.0)

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
    @test Base.isempty(GI.PhysicalGeometries.detect_intersection(
        triangle,
        SVector(0.25, 0.25, 1.0),
        SVector(0.25, 0.25, -1.0),
        triangle_hit,
    ))

    scale = Scale(geometry_3d, 2.0)
    @test scale isa Transformations.Transformation{3, 3, typeof(geometry_3d)}
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
    reflection = Rotate(geometry_2d, SMatrix{2, 2, Float64}([-1.0 0.0; 0.0 1.0]))
    @test Transformations.forward(reflection, position_2d) == SVector(-1.0, 2.0)

    projection = Rotate(geometry_2d, SMatrix{2, 3, Float64}([1.0 0.0 0.0; 0.0 1.0 0.0]))
    @test projection isa Transformations.Transformation{3, 2, typeof(geometry_2d)}
    @test Transformations.forward(projection, SVector(1.0, 2.0, 3.0)) == SVector(1.0, 2.0)
    @test_throws ArgumentError Transformations.backward(projection, SVector(1.0, 2.0))
    @test_throws ArgumentError Transformations.forward_normal(projection, SVector(1.0, 2.0, 3.0))
    @test Transformations.backward_normal(projection, SVector(1.0, 2.0)) ≈ SVector(1.0, 2.0, 0.0) / sqrt(5)
    @test_throws ArgumentError BoundingBoxes.InternalBoundingBox(projection)

    geometry_1d = TestGeometry{1}(0)
    projection_z = Rotate(geometry_1d, SMatrix{1, 3, Float64}([0.0 0.0 1.0]))
    @test Transformations.forward(projection_z, SVector(1.0, 2.0, 3.0)) == SVector(3.0)
    @test_throws ArgumentError Rotate(geometry_1d, SMatrix{1, 2, Float64}([2.0 0.0]))

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

    grid_box = BoundingBoxes.InternalBoundingBox([2.0, 2.0])
    grid = BoundingBoxes.grid_indices(
        grid_box,
        [4, 4],
        [
            BoundingBoxes.InternalBoundingBox([0.5, 0.5], [-1.5, -1.5]),
            BoundingBoxes.InternalBoundingBox([1.0, 0.5], [0.0, 0.0]),
        ],
    )
    @test grid[1, 1] == [1]
    @test grid[2, 2] == [2]
    @test grid[2, 3] == [2]
    @test grid[3, 2] == [2]
    @test grid[3, 3] == [2]
    @test grid[4, 4] == Int[]
    edge_grid = BoundingBoxes.grid_indices(
        BoundingBoxes.InternalBoundingBox([2.0]),
        [4],
        [
            BoundingBoxes.InternalBoundingBox([0.0], [-2.0]),
            BoundingBoxes.InternalBoundingBox([0.0], [-1.0]),
            BoundingBoxes.InternalBoundingBox([0.0], [2.0]),
        ],
    )
    @test edge_grid[1] == [1, 2]
    @test edge_grid[2] == [2]
    @test edge_grid[3] == Int[]
    @test edge_grid[4] == [3]
    @test_throws ArgumentError BoundingBoxes.grid_indices(grid_box, [0, 4], BoundingBoxes.InternalBoundingBox{2}[])
    @test_throws DimensionMismatch BoundingBoxes.grid_indices(grid_box, [4], BoundingBoxes.InternalBoundingBox{2}[])

    repeating_box = BoundingBoxes.InternalBoundingBox([3.0])
    shifts, repeating_grid = BoundingBoxes.grid_indices_repeating(
        BoundingBoxes.InternalBoundingBox([2.0]),
        [4],
        [4.0],
        [repeating_box],
    )
    @test shifts == [SVector(4.0), SVector(-4.0)]
    @test repeating_grid[1] == [(Int32(1), Int32(0)), (Int32(1), Int32(2))]
    @test repeating_grid[2] == [(Int32(1), Int32(0)), (Int32(1), Int32(2))]
    @test repeating_grid[3] == [(Int32(1), Int32(1)), (Int32(1), Int32(0))]
    @test repeating_grid[4] == [(Int32(1), Int32(1)), (Int32(1), Int32(0))]
    repeating_zero_boundary = BoundingBoxes.grid_indices_repeating(
        BoundingBoxes.InternalBoundingBox([2.0]),
        [4],
        [4.0],
        [BoundingBoxes.InternalBoundingBox([0.0], [-1.0])],
    )
    @test repeating_zero_boundary[1] == SVector{1, Float64}[]
    @test repeating_zero_boundary[2][1] == [(Int32(1), Int32(0))]
    @test repeating_zero_boundary[2][2] == [(Int32(1), Int32(0))]
    @test repeating_zero_boundary[2][3] == Tuple{Int32, Int32}[]
    @test repeating_zero_boundary[2][4] == Tuple{Int32, Int32}[]

    repeating_zero_lower_edge = BoundingBoxes.grid_indices_repeating(
        BoundingBoxes.InternalBoundingBox([2.0]),
        [4],
        [4.0],
        [BoundingBoxes.InternalBoundingBox([0.0], [-2.0])],
    )
    @test repeating_zero_lower_edge[2][1] == [(Int32(1), Int32(0))]
    @test all(isempty, repeating_zero_lower_edge[2][2:end])

    repeating_zero_upper_edge = BoundingBoxes.grid_indices_repeating(
        BoundingBoxes.InternalBoundingBox([2.0]),
        [4],
        [4.0],
        [BoundingBoxes.InternalBoundingBox([0.0], [2.0])],
    )
    @test repeating_zero_upper_edge[2][end] == [(Int32(1), Int32(0))]
    @test all(isempty, repeating_zero_upper_edge[2][1:end-1])
    @test_throws ArgumentError BoundingBoxes.grid_indices_repeating(
        BoundingBoxes.InternalBoundingBox([2.0]),
        [4],
        [0.0],
        BoundingBoxes.InternalBoundingBox{1}[],
    )

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

@testset "Fixed geometry susceptibility state" begin
    walls = mr.fix(mr.Walls(position=0.0))
    @test walls.susceptibility == ()
    susceptibility_position = SVector(0.0, 0.0, 0.0)
    @test GI.susceptibility_off_resonance(walls, susceptibility_position) == 0.0
    @test GI.susceptibility_off_resonance(walls, susceptibility_position, true) == 0.0
    @test GI.off_resonance_gradient(walls, 1.0) == 0.0

    cylinder = mr.Cylinders(
        radius=1.0,
        g_ratio=0.8,
        susceptibility_iso=1.0,
        susceptibility_aniso=0.0,
        grid_resolution=Inf,
    )
    fixed_cylinder = mr.fix(cylinder)
    @test fixed_cylinder.susceptibility isa Tuple
    @test length(fixed_cylinder.susceptibility) == 1
    @test fixed_cylinder.susceptibility[1] isa
        GI.Susceptibility.Grid.SusceptibilityGridNoRepeat
    @test GI.susceptibility_off_resonance(fixed_cylinder, susceptibility_position) ==
        GI.Susceptibility.susceptibility_off_resonance(
            fixed_cylinder.susceptibility,
            susceptibility_position,
        )
    @test GI.susceptibility_off_resonance(fixed_cylinder, susceptibility_position, true) ==
        GI.Susceptibility.susceptibility_off_resonance(
            fixed_cylinder.susceptibility,
            susceptibility_position,
            true,
        )
    @test GI.off_resonance_gradient(fixed_cylinder, 1.0) ==
        GI.Susceptibility.off_resonance_gradient(fixed_cylinder.susceptibility, 1.0)

    annulus = mr.Annuli(
        inner=0.7,
        outer=1.0,
        myelin=true,
        susceptibility_iso=1.0,
        susceptibility_aniso=0.0,
        repeats=[4.0, 4.0],
        grid_resolution=Inf,
    )
    fixed_annulus = mr.fix(annulus)
    @test fixed_annulus.susceptibility[1] isa
        GI.Susceptibility.Grid.SusceptibilityGridRepeat

    combined = mr.fix([cylinder, cylinder])
    @test length(combined.susceptibility) == 2
    @test all(
        susceptibility isa GI.Susceptibility.Grid.SusceptibilityGridNoRepeat
        for susceptibility in combined.susceptibility
    )
end
