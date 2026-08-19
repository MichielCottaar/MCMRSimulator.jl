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
const BaseObstructions = GI.PhysicalGeometries.BaseObstructions
const Meshes = GI.PhysicalGeometries.Meshes
const Mesh = Meshes.Mesh
const IntersectionGrid = GI.PhysicalGeometries.GridDispatch.IntersectionGrid
const Transparents = GI.PhysicalGeometries.Transparents
const Transparent = GI.PhysicalGeometries.Transparent
const SizeScaleOverride = GI.SizeScaleOverride
const size_scale = GI.size_scale
const has_inside = GI.PhysicalGeometries.has_inside
const has_single_inside = GI.PhysicalGeometries.has_single_inside
const isinside_single = GI.PhysicalGeometries.isinside_single
const inside_indices = GI.PhysicalGeometries.inside_indices
const Properties = GI.Properties
const Plot = mr.Plot

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

@testset "intersection grids" begin
    grid_box = BoundingBoxes.InternalBoundingBox([2.0, 1.0, 1.0])
    child_boxes = [
        BoundingBoxes.InternalBoundingBox([0.5, 1.0, 1.0], [-1.0, 0.0, 0.0]),
        BoundingBoxes.InternalBoundingBox([0.5, 1.0, 1.0], [1.0, 0.0, 0.0]),
    ]
    grid = IntersectionGrid(child_boxes; resolution=2.0, bounding_box=grid_box)
    @test grid.bounding_box == grid_box
    @test size(grid.indices) == (2, 1, 1)
    @test grid.indices[1, 1, 1] == [1]
    @test grid.indices[2, 1, 1] == [2]

    planar_box = BoundingBoxes.InternalBoundingBox([1.0, 1.0, 0.0])
    planar_grid = IntersectionGrid([planar_box]; resolution=1.0)
    @test planar_grid.bounding_box == planar_box
    @test all(isfinite, planar_grid.inv_resolution)
    @test BoundingBoxes.upper(planar_grid.grid_bounding_box)[3] > BoundingBoxes.lower(planar_grid.grid_bounding_box)[3]
end

@testset "projected mesh field of view" begin
    plot_plane = mr.PlotPlane(sizex=1.0, sizey=1.0)
    crossing_mesh = [(;
        vertices=[
            SVector(-2.0, 0.0, -1.0),
            SVector(2.0, 0.0, -1.0),
            SVector(0.0, 0.0, 1.0),
        ],
        triangles=[(1, 2, 3)],
    )]
    projected = Plot.project_mesh_plane(plot_plane, crossing_mesh)
    @test projected[1:2] == [SVector(0.5, 0.0), SVector(-0.5, 0.0)]
    @test all(isnan, projected[3])

    outside_mesh = [(;
        vertices=[
            SVector(-2.0, 2.0, -1.0),
            SVector(2.0, 2.0, -1.0),
            SVector(0.0, 2.0, 1.0),
        ],
        triangles=[(1, 2, 3)],
    )]
    @test isempty(Plot.project_mesh_plane(plot_plane, outside_mesh))

    coplanar_mesh = [(;
        vertices=[
            SVector(-2.0, 0.0, 0.0),
            SVector(2.0, 0.0, 0.0),
            SVector(0.0, 2.0, 0.0),
        ],
        triangles=[(1, 2, 3)],
    )]
    coplanar_projected = Plot.project_mesh_plane(plot_plane, coplanar_mesh)
    finite_points = filter(point -> all(isfinite, point), coplanar_projected)
    @test all(all(-0.5 <= coordinate <= 0.5 for coordinate in point) for point in finite_points)
    @test any(point -> point == SVector(-0.5, 0.0), finite_points)
    @test any(point -> point == SVector(0.5, 0.0), finite_points)
end

const get_value = Properties.get_value
const all_property_values = Properties.all_property_values

@testset "reflections" begin
    Reflections = mr.Reflections
    index = (1, 2)
    collision = GI.Intersection(0.5, index, index, SVector(1.0, 0.0, 0.0), false, false)

    free = Reflections.Reflection(1.0)
    @test !Reflections.has_intersection(free)
    @test isnothing(Reflections.previous_hit(free))

    reflection = Reflections.Reflection(
        collision,
        SVector(-1.0, 0.0, 0.0),
        1.0,
        0.0,
        0.0,
    )
    @test Reflections.has_intersection(reflection)
    @test Reflections.has_hit(reflection) == index
    @test Reflections.previous_hit(reflection) === collision
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
    @test empty_fixed.volume.R1.value == 0.0
    @test empty_fixed.volume.R2.value == 0.0
    @test empty_fixed.volume.off_resonance.value == 0.0
    @test empty_fixed.surface.permeability.value == 0.0
    @test empty_fixed.surface.dwell_time.value == 0.0
    @test empty_fixed.surface.surface_relaxation.value == 0.0
    @test empty_fixed.surface.density.value == 0.0
    configured_empty = mr.fix([]; permeability=1.0, dwell_time=2.0, surface_relaxation=3.0, density=4.0)
    @test configured_empty.surface.permeability.value == 1.0
    @test configured_empty.surface.dwell_time.value == 2.0
    @test configured_empty.surface.surface_relaxation.value == 3.0
    @test configured_empty.surface.density.value == 4.0

    property_geometry = mr.fix(
        mr.Spheres(
            radius=[1.0, 2.0],
            permeability_surface=[3.0, 4.0],
            density_surface=[5.0, 6.0],
            dwell_time_surface=[7.0, 8.0],
            relaxation_surface=[9.0, 10.0],
        ),
    )
    property_intersection = GI.Intersection(
        0.5,
        (2, 1),
        (2, 1),
        SVector(1.0, 0.0, 0.0),
        false,
        false,
    )
    @test GI.permeability(property_geometry, property_intersection) == 4.0
    @test GI.surface_relaxation(property_geometry, property_intersection) == 10.0
    @test GI.surface_density(property_geometry, property_intersection) == 6.0
    @test GI.dwell_time(property_geometry, property_intersection) == 8.0
    gap_intersection = GI.Intersection(
        property_intersection.distance,
        property_intersection.indices,
        property_intersection.property_indices,
        property_intersection.normal,
        property_intersection.inside,
        true,
    )
    @test GI.permeability(property_geometry, gap_intersection) == Inf
    @test GI.surface_relaxation(property_geometry, gap_intersection) == 0.0
    @test GI.surface_density(property_geometry, gap_intersection) == 0.0
    @test GI.dwell_time(property_geometry, gap_intersection) == 0.0
    mri_geometry = mr.fix(mr.Spheres(
        radius=[1.0, 2.0],
        R1_inside=[1.0, 2.0],
        R2_inside=[3.0, 4.0],
        off_resonance_inside=[5.0, 6.0],
        R1_surface=[10.0, 20.0],
        R2_surface=[30.0, 40.0],
        off_resonance_surface=[50.0, 60.0],
    ))
    inside_second = GI.IsInside(SVector{1, Tuple{Int}}(((2,),)))
    mri_intersection = GI.Intersection(
        0.5,
        (2,),
        (2,),
        SVector(1.0, 0.0, 0.0),
        false,
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
    inside_both = GI.IsInside(SVector{2, Tuple{Int}}((1,), (2,)))
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
    equal_depth_tuple = GeometryTuple{3}((
        BaseObstructions.Sphere(1.0),
        BaseObstructions.Sphere(2.0),
    ))
    equal_tuple_inside = @inferred inside_indices(equal_depth_tuple, SVector(0.0, 0.0, 0.0))
    @test equal_tuple_inside isa Vector{Tuple{Int}}
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
    @test shifted_render_mesh[1].vertices == [vertex + SVector(1.0, 2.0, 3.0) for vertex in mesh.vertices]
    grouped_render_mesh = GI.geometry_mesh(GeometryTuple{3}((mesh, mesh)))
    @test length(grouped_render_mesh) == 2

    wall = BaseObstructions.InfiniteWall()
    wall_rotation = Rotate(wall, reshape([1.0, 0.0, 0.0], 3, 1))
    wall_mesh = GI.geometry_mesh(wall_rotation; height=2.0)
    @test length(wall_mesh) == 1
    @test length(wall_mesh[1].vertices) == 4
    @test all(vertex[1] == 0.0 for vertex in wall_mesh[1].vertices)
    wall_bounded_mesh = GI.geometry_mesh(
        wall_rotation,
        PublicBoundingBoxes.BoundingBox([-1.0, -2.0, -3.0], [1.0, 2.0, 3.0]),
    )
    wall_extents = [
        maximum(vertex[dimension] for vertex in wall_bounded_mesh[1].vertices) -
        minimum(vertex[dimension] for vertex in wall_bounded_mesh[1].vertices)
        for dimension in 2:3
    ]
    @test sort(wall_extents) == [4.0, 6.0]

    repeated_walls = Repeat(GeometryVector{1}([wall]), [4.0])
    repeated_wall_mesh = GI.geometry_mesh(
        Rotate(repeated_walls, reshape([1.0, 0.0, 0.0], 3, 1)),
        PublicBoundingBoxes.BoundingBox([-5.0, -2.0, -2.0], [5.0, 2.0, 2.0]);
        height=2.0,
    )
    @test length(repeated_wall_mesh) == 3

    cylinder = BaseObstructions.InfiniteCylinder(1.0)
    cylinder_rotation = Rotate(cylinder, [1.0 0.0; 0.0 1.0; 0.0 0.0])
    cylinder_mesh = GI.geometry_mesh(cylinder_rotation; nsamples=8, height=2.0)
    @test length(cylinder_mesh) == 1
    @test length(cylinder_mesh[1].vertices) == 16
    @test length(cylinder_mesh[1].triangles) == 16
    cylinder_bounded_mesh = GI.geometry_mesh(
        cylinder_rotation,
        PublicBoundingBoxes.BoundingBox([-1.0, -1.0, -3.0], [1.0, 1.0, 3.0]);
        nsamples=8,
    )
    @test maximum(vertex[3] for vertex in cylinder_bounded_mesh[1].vertices) -
        minimum(vertex[3] for vertex in cylinder_bounded_mesh[1].vertices) == 6.0

    sphere_mesh = GI.geometry_mesh(BaseObstructions.Sphere(1.0); nsamples=20)
    @test length(sphere_mesh) == 1
    @test length(sphere_mesh[1].triangles) == 36

    repeated_circles = Repeat(GeometryVector{2}([cylinder]), [4.0, 4.0])
    repeated_cylinder_mesh = GI.geometry_mesh(
        Rotate(repeated_circles, [1.0 0.0; 0.0 1.0; 0.0 0.0]),
        PublicBoundingBoxes.BoundingBox([-5.0, -5.0, -1.0], [5.0, 5.0, 1.0]);
        nsamples=8,
        height=2.0,
    )
    @test length(repeated_cylinder_mesh) == 9

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
    @test has_single_inside(typeof(mesh))
    @test isinside_single(mesh, SVector(0.2, 0.2, 0.2))
    @test !isinside_single(mesh, SVector(1.2, 0.2, 0.2))
    stuck_inside_mesh = (1, true, 0.5)
    stuck_outside_mesh = (1, false, 0.5)
    @test isinside_single(mesh, SVector(1.2, 0.2, 0.2), stuck_inside_mesh)
    @test !isinside_single(mesh, SVector(0.2, 0.2, 0.2), stuck_outside_mesh)
    @test mesh.indices[mesh.first_index_of_gap] == SVector(4, 2, 3)
    @test BoundingBoxes.lower(mesh.bounding_box) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.upper(mesh.bounding_box) == SVector(1.0, 1.0, 1.0)
    stability_start = SVector(-2.0, 0.0, 0.0)
    stability_destination = SVector(2.0, 0.0, 0.0)
    stable_group = GeometryVector([BaseObstructions.Sphere(1.0)])
    stable_group_hit = GI.PhysicalGeometries.find_intersection(
        stable_group, stability_start, stability_destination,
    )
    @test stable_group_hit isa Tuple
    stable_group_inside = @inferred inside_indices(stable_group, SVector(0.0, 0.0, 0.0))
    @test stable_group_inside isa Vector{Tuple{Int}}
    stable_mesh_hit = GI.PhysicalGeometries.find_intersection(
        mesh, stability_start, stability_destination,
    )
    @test stable_mesh_hit isa Tuple
    stable_mesh_group = GeometryVector([mesh])
    stable_mesh_group_hit = GI.PhysicalGeometries.find_intersection(
        stable_mesh_group, stability_start, stability_destination,
    )
    @test stable_mesh_group_hit isa Tuple
    stable_mesh_group_inside = @inferred inside_indices(stable_mesh_group, SVector(0.2, 0.2, 0.2))
    @test stable_mesh_group_inside isa Vector{Tuple{Int}}
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
    mesh_real_hit = GI.PhysicalGeometries.find_intersection(
        mesh,
        SVector(0.2, 0.2, -1.0),
        SVector(0.2, 0.2, 0.2),
    )
    @test mesh_real_hit[3] ≈ 5 / 6
    @test mesh_real_hit[1] == 1
    mesh_gap_hit = GI.PhysicalGeometries.find_intersection(
        mesh,
        SVector(0.2, 0.2, 0.2),
        SVector(1.0, 1.0, 1.0),
    )
    @test mesh_gap_hit[3] ≈ 1 / 6
    @test mesh_gap_hit[1] == 4
    @test isnothing(GI.PhysicalGeometries.find_intersection(
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
    @test inside_indices(repeated_sphere, SVector(4.5, 0.0, 0.0)) == [(SVector(1, 0, 0),)]
    @test inside_indices(repeated_sphere, SVector(2.0, 0.0, 0.0)) == Tuple{SVector{3, Int}}[]
    shifted_repeated_sphere = Repeat(
        Shift(BaseObstructions.Sphere(0.6), SVector(-1.0, 0.0, 0.0)),
        [2.0, 4.0, 4.0],
    )
    shifted_position = SVector(-0.5, 0.0, 0.0)
    @test inside_indices(shifted_repeated_sphere, shifted_position) == [(SVector(0, 0, 0),)]
    neighboring_position = SVector(0.5, 0.0, 0.0)
    @test Repeats._candidate_shifts(shifted_repeated_sphere, neighboring_position) == [
        SVector(0, 0, 0), SVector(-1, 0, 0),
    ]
    @test inside_indices(shifted_repeated_sphere, neighboring_position) == [(SVector(1, 0, 0),)]
    @test_throws ArgumentError Repeat(BaseObstructions.Sphere(1.0), [0.0, 4.0, 4.0])
    @test_throws ArgumentError Repeat(BaseObstructions.Sphere(2.0), [1.0, 4.0, 4.0])
    @test_throws ArgumentError Repeat(
        Shift(BaseObstructions.Sphere(1.0), [3.0, 0.0, 0.0]),
        [2.0, 4.0, 4.0],
    )
    @test_throws ArgumentError BoundingBoxes.InternalBoundingBox(repeated_sphere)
    repeated_hit = GI.PhysicalGeometries.find_intersection(
        repeated_sphere,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    @test repeated_hit[end] ≈ 0.4
    overlapping_repeat = Repeat(BaseObstructions.Sphere(1.5), [2.0, 4.0, 4.0])
    @test overlapping_repeat.lower_overlap == SVector(0.5, 0.0, 0.0)
    @test overlapping_repeat.upper_overlap == SVector(0.5, 0.0, 0.0)
    @test inside_indices(overlapping_repeat, SVector(1.4, 0.0, 0.0)) == [
        (SVector(1, 0, 0),), (SVector(0, 0, 0),),
    ]
    overlapping_hit = GI.PhysicalGeometries.find_intersection(
        overlapping_repeat,
        SVector(1.4, 0.0, 0.0),
        SVector(1.9, 0.0, 0.0),
    )
    @test overlapping_hit[end] ≈ 0.2

    @test get_value(3.0, (1, 2)) == 3.0
    leaf_properties = Properties.GeometryLeafProperties(3.0)
    @test leaf_properties isa Properties.GeometryProperties{Float64}
    @test eltype(leaf_properties) === Float64
    vector_properties = Properties.GeometryVectorProperties([10.0, 20.0])
    @test vector_properties isa Properties.GeometryProperties{Float64}
    @test eltype(vector_properties) === Float64
    @test get_value(vector_properties, (2, 1)) == 20.0
    tuple_properties = Properties.GeometryTupleProperties((1.0, vector_properties))
    @test tuple_properties isa Properties.GeometryProperties{Float64}
    @test eltype(tuple_properties) === Float64
    @test get_value(tuple_properties, (2, 2)) == 20.0
    @test get_value(tuple_properties, [
        (1, 4),
        (2, 2),
    ]) == 21.0
    aggregation_properties = Properties.GeometryTupleProperties((
        Properties.GeometryVectorProperties([1.0, 2.0]),
        10.0,
    ))
    @test get_value(aggregation_properties, [
        (1, 1),
        (1, 2),
        (2, 1),
        (2, 2),
    ]) == 13.0
    @test get_value(aggregation_properties, Tuple[]) == 0
    @test_throws ArgumentError get_value(aggregation_properties, [
        (2, 1),
        (1, 1),
    ])
    @test_throws BoundsError get_value(vector_properties, (3,))
    @test_throws ArgumentError Properties.GeometryVectorProperties([
        Properties.GeometryLeafProperties(1.0),
        Properties.GeometryLeafProperties(2),
    ])
    @test_throws ArgumentError Properties.GeometryTupleProperties((
        Properties.GeometryLeafProperties(1.0),
        Properties.GeometryLeafProperties(2),
    ))

    sphere = BaseObstructions.Sphere(1.0)
    @test has_single_inside(typeof(sphere))
    @test isinside_single(sphere, SVector(0.0, 0.0, 0.0))
    @test !isinside_single(sphere, SVector(2.0, 0.0, 0.0))
    @test !has_inside(typeof(BaseObstructions.InfiniteWall()))

    shifted_inside = Shift(sphere, [2.0, 0.0, 0.0])
    @test has_single_inside(typeof(shifted_inside))
    @test isinside_single(shifted_inside, SVector(2.0, 0.0, 0.0))
    @test !isinside_single(shifted_inside, SVector(0.0, 0.0, 0.0))

    inside_spheres = GeometryVector{3}([
        Shift(sphere, [-0.5, 0.0, 0.0]),
        Shift(sphere, [0.5, 0.0, 0.0]),
    ])
    @test inside_indices(inside_spheres, SVector(0.0, 0.0, 0.0)) == [
        (1,),
        (2,),
    ]
    walls = GeometryVector{1}([BaseObstructions.InfiniteWall(), BaseObstructions.InfiniteWall()])
    @test inside_indices(walls, SVector(0.0)) == Tuple{Int}[]

    nested = GeometryTuple{3}((inside_spheres, sphere))
    @test inside_indices(nested, SVector(0.0, 0.0, 0.0)) == [
        (1, 1),
        (1, 2),
        (2,),
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
    @test Transformations.to_child_coordinates(shift, position) == SVector(3.0, 3.0, 3.0)
    @test Transformations.from_child_coordinates(shift, SVector(3.0, 3.0, 3.0)) == position
    @test Transformations.to_child_coordinates_normal(shift, normal) == normal
    @test Transformations.from_child_coordinates_normal(shift, normal) == normal

    infinite_wall = BaseObstructions.InfiniteWall()
    @test infinite_wall isa PhysicalGeometry{1}
    @test BoundingBoxes.InternalBoundingBox(infinite_wall) isa BoundingBoxes.InternalBoundingBox{1}
    wall_hit = GI.PhysicalGeometries.find_intersection(
        infinite_wall,
        SVector(-1.0),
        SVector(1.0),
    )
    @test wall_hit[2] == 0.5
    @test GI.PhysicalGeometries.get_intersection_params(
        infinite_wall, SVector(-1.0), SVector(1.0), wall_hit,
    ).normal == SVector(-1.0)
    previous_wall_hit = (false, 0.5)
    @test isnothing(GI.PhysicalGeometries.find_intersection(
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
    sphere_hit = GI.PhysicalGeometries.find_intersection(
        sphere,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
    )
    @test sphere_hit[2] ≈ 1 / 6
    @test GI.PhysicalGeometries.get_intersection_params(
        sphere, SVector(-3.0, 0.0, 0.0), SVector(3.0, 0.0, 0.0), sphere_hit,
    ).normal == SVector(-1.0, 0.0, 0.0)
    boundary = SVector(2.0, 0.0, 0.0)
    boundary_inside = (true, 0.5)
    boundary_outside = (false, 0.5)
    @test isinside_single(sphere, boundary, boundary_inside)
    @test !isinside_single(sphere, boundary, boundary_outside)

    fixed_sphere = mr.fix(mr.Spheres(radius=2.0))
    fixed_boundary_inside = (1, true, 0.5)
    fixed_boundary_outside = (1, false, 0.5)
    @test GI.isinside(fixed_sphere, boundary, fixed_boundary_inside).inside_of == [(1,)]
    @test GI.isinside(fixed_sphere, boundary, fixed_boundary_outside).inside_of == Tuple{Int}[]

    @test isnothing(GI.PhysicalGeometries.find_intersection(
        sphere,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
        sphere_hit,
    ))
    previous_inside_sphere = (true, 0.5)
    inside_sphere_hit = GI.PhysicalGeometries.find_intersection(
        sphere,
        SVector(-1.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
        previous_inside_sphere,
    )
    @test inside_sphere_hit[2] ≈ 3 / 4
    @test GI.PhysicalGeometries.get_intersection_params(
        sphere, SVector(-1.0, 0.0, 0.0), SVector(3.0, 0.0, 0.0), inside_sphere_hit,
    ).normal == SVector(-1.0, 0.0, 0.0)

    fixed_sphere_hit = GI.detect_intersection(
        fixed_sphere,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
    )
    underlying_sphere_hit = GI.PhysicalGeometries.find_intersection(
        fixed_sphere.geometry,
        SVector(-3.0, 0.0, 0.0),
        SVector(3.0, 0.0, 0.0),
    )
    @test fixed_sphere_hit.distance == underlying_sphere_hit[end]

    shifted_sphere = Shift(sphere, [1.0, 0.0, 0.0])
    shifted_sphere_box = BoundingBoxes.InternalBoundingBox(shifted_sphere)
    @test BoundingBoxes.lower(shifted_sphere_box) == SVector(-1.0, -2.0, -2.0)
    @test BoundingBoxes.upper(shifted_sphere_box) == SVector(3.0, 2.0, 2.0)
    shifted_hit = GI.PhysicalGeometries.find_intersection(
        shifted_sphere,
        SVector(-4.0, 0.0, 0.0),
        SVector(2.0, 0.0, 0.0),
    )
    @test shifted_hit[end] ≈ 1 / 2
    @test GI.PhysicalGeometries.get_intersection_params(
        shifted_sphere, SVector(-4.0, 0.0, 0.0), SVector(2.0, 0.0, 0.0), shifted_hit,
    ).normal == SVector(-1.0, 0.0, 0.0)

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
    filtered_hit = GI.PhysicalGeometries.find_intersection(
        filtered_spheres,
        SVector(-2.0, 0.0, 0.0),
        SVector(2.0, 0.0, 0.0),
    )
    @test filtered_hit[1] == 1
    @test inside_indices(filtered_spheres, SVector(0.0, 0.0, 0.0)) == [
        (1,),
    ]
    grid_spheres = GeometryVector([
        Shift(sphere, [-3.0, 0.0, 0.0]),
        Shift(sphere, [3.0, 0.0, 0.0]),
    ]; grid=true, grid_resolution=2.0)
    @test grid_spheres isa GeometryVectorGrid{3, Shift{3, BaseObstructions.Sphere}}
    @test inside_indices(grid_spheres, SVector(3.0, 0.0, 0.0)) == [(2,)]
    regular_grid_hit = GI.PhysicalGeometries.find_intersection(
        GeometryVector([
            Shift(sphere, [-3.0, 0.0, 0.0]),
            Shift(sphere, [3.0, 0.0, 0.0]),
        ]),
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    grid_hit = GI.PhysicalGeometries.find_intersection(
        grid_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    @test grid_hit == regular_grid_hit
    @test_throws ArgumentError GeometryVectorBoundingBox(
        [sphere],
        BoundingBoxes.InternalBoundingBox{3}[],
    )
    @test_throws ArgumentError BoundingBoxes.InternalBoundingBox(GeometryVector{3}(TestGeometry{3}[]))
    @test_throws ArgumentError BoundingBoxes.InternalBoundingBox(GeometryTuple{3}(()))
    grouped_hit = GI.PhysicalGeometries.find_intersection(
        grouped_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    @test grouped_hit[end] ≈ 1 / 5
    @test grouped_hit[1] == 1
    grouped_boundary = SVector(3.0, 0.0, 0.0)
    grouped_inside = (2, true, 0.5)
    grouped_outside = (2, false, 0.5)
    @test inside_indices(grouped_spheres, grouped_boundary, grouped_inside) == [(2,)]
    @test inside_indices(grouped_spheres, grouped_boundary, grouped_outside) == Tuple{Int}[]

    next_grouped_hit = GI.PhysicalGeometries.find_intersection(
        grouped_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
        grouped_hit,
    )
    @test next_grouped_hit[end] ≈ 3 / 5
    @test next_grouped_hit[1] == 2
    @test GI.PhysicalGeometries.find_intersection(
        grouped_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
        (0, false, 0.5),
    ) == grouped_hit
    @test GI.PhysicalGeometries.find_intersection(
        grouped_spheres,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
        (3, false, 0.5),
    ) == grouped_hit

    scaled_sphere = Scale(sphere, 2.0)
    scaled_hit = GI.PhysicalGeometries.find_intersection(
        scaled_sphere,
        SVector(-6.0, 0.0, 0.0),
        SVector(6.0, 0.0, 0.0),
    )
    @test scaled_hit[end] ≈ 1 / 6
    @test scaled_hit[1] == false

    projected_cylinder = Rotate(cylinder, SMatrix{3, 2, Float64}([1.0 0.0; 0.0 1.0; 0.0 0.0]))
    projected_hit = GI.PhysicalGeometries.find_intersection(
        projected_cylinder,
        SVector(-3.0, 0.0, 4.0),
        SVector(3.0, 0.0, 4.0),
    )
    @test projected_hit[end] ≈ 1 / 6
    @test projected_hit[1] == false

    triangle = BaseObstructions.FullTriangle(
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
    )
    triangle_box = BoundingBoxes.InternalBoundingBox(triangle)
    @test triangle_box isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.lower(triangle_box) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.upper(triangle_box) == SVector(1.0, 1.0, 0.0)
    triangle_hit = GI.PhysicalGeometries.find_intersection(
        triangle,
        SVector(0.25, 0.25, 1.0),
        SVector(0.25, 0.25, -1.0),
    )
    @test triangle_hit[end] == 0.5
    @test triangle_hit[1] == false
    @test isnothing(GI.PhysicalGeometries.find_intersection(
        triangle,
        SVector(0.25, 0.25, 1.0),
        SVector(0.25, 0.25, -1.0),
        triangle_hit,
    ))

    scale = Scale(geometry_3d, 2.0)
    @test scale isa Transformations.Transformation{3, 3, typeof(geometry_3d)}
    @test scale.scale == 2.0
    @test Transformations.to_child_coordinates(scale, position) == SVector(2.0, 2.5, 3.0)
    @test Transformations.from_child_coordinates(scale, SVector(2.0, 2.5, 3.0)) == position
    @test GI.PhysicalGeometries.to_child_coordinates_normal(scale, normal) == normal
    @test GI.PhysicalGeometries.from_child_coordinates_normal(scale, normal) == normal
    @test_throws ArgumentError Scale(geometry_3d, 0.0)
    @test_throws ArgumentError Scale(geometry_3d, -1.0)
    rotation = Rotate(geometry_2d, SMatrix{2, 2, Float64}([0.0 -1.0; 1.0 0.0]))
    position_2d = SVector(1.0, 2.0)
    @test Transformations.to_child_coordinates(rotation, position_2d) == SVector(2.0, -1.0)
    @test Transformations.from_child_coordinates(rotation, SVector(2.0, -1.0)) == position_2d
    @test GI.PhysicalGeometries.to_child_coordinates_normal(rotation, position_2d) == SVector(2.0, -1.0)
    reflection = Rotate(geometry_2d, SMatrix{2, 2, Float64}([-1.0 0.0; 0.0 1.0]))
    @test Transformations.to_child_coordinates(reflection, position_2d) == SVector(-1.0, 2.0)

    projection = Rotate(geometry_2d, SMatrix{3, 2, Float64}([1.0 0.0; 0.0 1.0; 0.0 0.0]))
    @test projection isa Transformations.Transformation{3, 2, typeof(geometry_2d)}
    @test Transformations.to_child_coordinates(projection, SVector(1.0, 2.0, 3.0)) == SVector(1.0, 2.0)
    @test Transformations.from_child_coordinates(projection, SVector(1.0, 2.0)) == SVector(1.0, 2.0, 0.0)
    @test_throws ArgumentError GI.PhysicalGeometries.to_child_coordinates_normal(projection, SVector(1.0, 2.0, 3.0))
    @test GI.PhysicalGeometries.from_child_coordinates_normal(projection, SVector(1.0, 2.0)) ≈ SVector(1.0, 2.0, 0.0) / sqrt(5)
    @test_throws ArgumentError BoundingBoxes.InternalBoundingBox(projection)

    geometry_1d = TestGeometry{1}(0)
    projection_z = Rotate(geometry_1d, SMatrix{3, 1, Float64}([0.0; 0.0; 1.0]))
    @test Transformations.to_child_coordinates(projection_z, SVector(1.0, 2.0, 3.0)) == SVector(3.0)
    @test_throws ArgumentError Rotate(geometry_1d, SMatrix{2, 1, Float64}([2.0; 0.0]))

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
    scaled_centered_cube = Transformations.to_child_coordinates(scale, centered_cube)
    scaled_rect = Transformations.to_child_coordinates(scale, rect)
    @test scaled_centered_cube isa BoundingBoxes.InternalBoundingBox{3}
    @test scaled_rect isa BoundingBoxes.InternalBoundingBox{3}
    @test BoundingBoxes.lower(scaled_centered_cube) == SVector(-1.0, -1.0, -1.0)
    @test BoundingBoxes.upper(scaled_centered_cube) == SVector(1.0, 1.0, 1.0)
    @test BoundingBoxes.lower(scaled_rect) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.upper(scaled_rect) == SVector(1.0, 2.0, 3.0)
    @test Transformations.to_child_coordinates(shift, centered_cube) isa BoundingBoxes.InternalBoundingBox{3}
    transformed_centered_cube = Transformations.to_child_coordinates(shift, centered_cube)
    recovered_centered_cube = Transformations.from_child_coordinates(shift, transformed_centered_cube)
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
    @test Transformations.to_child_coordinates(projection, centered_cube) isa BoundingBoxes.InternalBoundingBox{2}
    @test Transformations.to_child_coordinates(projection_z, centered_cube) isa BoundingBoxes.InternalBoundingBox{1}
    @test_throws ArgumentError Transformations.from_child_coordinates(projection, centered_cube)

    rotation_3d = Rotate(geometry_3d, SMatrix{3, 3, Float64}(I))
    rotated_box = Transformations.to_child_coordinates(rotation_3d, centered_cube)
    @test rotated_box isa BoundingBoxes.InternalBoundingBox{3}
    transformed_bounds = (BoundingBoxes.lower(centered_cube), BoundingBoxes.upper(centered_cube))
    @test all(BoundingBoxes.isinside(rotated_box, Transformations.to_child_coordinates(rotation_3d, point)) for point in transformed_bounds)
    @test Transformations.from_child_coordinates(rotation_3d, rotated_box) isa BoundingBoxes.InternalBoundingBox{3}
end

@testset "Fixed geometry susceptibility state" begin
    empty_geometry = mr.fix(())
    @test GI.detect_intersection(empty_geometry, SVector(0.0, 0.0, 0.0), SVector(1.0, 0.0, 0.0)) === nothing
    @test length(GI.isinside(empty_geometry, SVector(0.0, 0.0, 0.0))) == 0
    @test GI.susceptibility_off_resonance(empty_geometry, SVector(0.0, 0.0, 0.0)) == 0.0

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
