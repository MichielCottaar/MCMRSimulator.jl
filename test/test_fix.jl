const FixBaseGeometry = mr.Geometries.Fix.FixBaseGeometry
const FixTransformations = mr.Geometries.Fix.FixTransformations
const FixProperties = mr.Geometries.Fix.FixProperties
const BaseObstructions = mr.Geometries.Internal.PhysicalGeometries.BaseObstructions
const Groups = mr.Geometries.Internal.PhysicalGeometries.Groups
const Meshes = mr.Geometries.Internal.PhysicalGeometries.Meshes
const Repeats = mr.Geometries.Internal.PhysicalGeometries.Repeats
const Transparent = mr.Geometries.Internal.PhysicalGeometries.Transparent
const SizeScaleOverride = mr.Geometries.Internal.SizeScaleOverride
const Transformations = mr.Geometries.Internal.PhysicalGeometries.Transformations
const Properties = mr.Geometries.Internal.Properties
const ObstructionIndex = mr.Geometries.Internal.Indices.ObstructionIndex
const size_scale = mr.Geometries.Internal.size_scale
const get_value = Properties.get_value
const fix = mr.Geometries.Fix.fix

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

    mesh = mr.Mesh(
        vertices=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
            [3.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
            [3.0, 1.0, 0.0],
            [3.0, 0.0, 1.0],
        ],
        triangles=[
            [1, 2, 3], [1, 2, 4], [1, 3, 4],
            [5, 6, 7], [5, 6, 8], [5, 7, 8],
        ],
    )
    mesh_geometry = FixBaseGeometry.fix_base_geometry(mesh)
    @test length(mesh_geometry) == 2
    @test all(mesh_part isa Meshes.Mesh for mesh_part in mesh_geometry)
    @test all(length(mesh_part.indices) == 4 for mesh_part in mesh_geometry)

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

    cylinders = mr.Cylinders(
        radius=[1.0, 2.0],
        R1_inside=[1.0, 2.0],
        R2_inside=3.0,
        off_resonance_inside=4.0,
        R1_surface=5.0,
        permeability_surface=[6.0, 7.0],
    )
    cylinder_properties = FixProperties.fix_properties(
        cylinders;
        dwell_time=8.0,
        surface_relaxation=9.0,
        density=10.0,
    )
    @test cylinder_properties.volume.R1 isa Properties.GeometryVectorProperties
    @test get_value(cylinder_properties.volume.R1, ObstructionIndex(SVector(2))) == 2.0
    @test get_value(cylinder_properties.volume.R2, ObstructionIndex(SVector(1))) == 3.0
    @test get_value(cylinder_properties.surface.permeability, ObstructionIndex(SVector(2))) == 7.0
    @test get_value(cylinder_properties.surface.dwell_time, ObstructionIndex(SVector(1))) == 8.0
    @test get_value(cylinder_properties.surface.surface_relaxation, ObstructionIndex(SVector(1))) == 9.0
    @test get_value(cylinder_properties.surface.density, ObstructionIndex(SVector(1))) == 10.0

    wall_properties = FixProperties.fix_properties(mr.Walls(position=0.0, R1_surface=2.0))
    @test get_value(wall_properties.volume.R1, ObstructionIndex(SVector(1))) == 0.0
    @test get_value(wall_properties.surface.R1, ObstructionIndex(SVector(1))) == 2.0

    annulus_properties = FixProperties.fix_properties(mr.Annuli(
        inner=1.0,
        outer=2.0,
        R1_inner_volume=1.0,
        R1_outer_volume=2.0,
        permeability_inner_surface=3.0,
        permeability_outer_surface=4.0,
    ))
    @test annulus_properties.volume.R1 isa Properties.GeometryTupleProperties
    @test get_value(annulus_properties.volume.R1, ObstructionIndex(SVector(1, 1))) == 1.0
    @test get_value(annulus_properties.volume.R1, ObstructionIndex(SVector(2, 1))) == 2.0
    @test get_value(annulus_properties.surface.permeability, ObstructionIndex(SVector(1, 1))) == 3.0
    @test get_value(annulus_properties.surface.permeability, ObstructionIndex(SVector(2, 1))) == 4.0
end

@testset "Fix dispatcher" begin
    group = mr.Spheres(
        radius=1.0,
        R1_inside=2.0,
        permeability_surface=3.0,
    )
    fixed = fix(group; dwell_time=4.0, surface_relaxation=5.0, density=6.0)
    @test fixed isa mr.Geometries.Fix.FixedGeometry
    @test fixed.volume.R1.value == 2.0
    @test fixed.surface.permeability.value == 3.0
    @test fixed.surface.dwell_time.value == 4.0
    @test fixed.surface.surface_relaxation.value == 5.0
    @test fixed.surface.density.value == 6.0
    @test fix(fixed) === fixed

    combined = fix([mr.Walls(position=0.0), group])
    @test combined.geometry isa mr.Geometries.Internal.PhysicalGeometries.Groups.GeometryTuple
    @test combined.volume.R1 isa Properties.GeometryTupleProperties
    @test combined.surface.density isa Properties.GeometryTupleProperties

    @test fix([]).geometry isa mr.Geometries.Internal.PhysicalGeometries.Groups.GeometryTuple
    @test_throws MethodError fix(group; R1=1.0)
    @test_throws MethodError fix(group; R2=1.0)
    @test_throws MethodError fix(group; off_resonance=1.0)
    @test_throws MethodError fix(group; relaxation=1.0)

    mesh = mr.Mesh(
        vertices=[
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
        triangles=[[1, 2, 3], [1, 2, 4], [1, 3, 4]],
        R1_inside=2.0,
        permeability_surface=3.0,
    )
    fixed_mesh = fix(mesh; dwell_time=4.0, surface_relaxation=5.0, density=6.0)
    @test fixed_mesh.geometry isa Groups.GeometryVectorBoundingBox
    @test first(fixed_mesh.geometry.geometries) isa Meshes.Mesh
    @test fixed_mesh.volume.R1.value == 2.0
    @test fixed_mesh.surface.permeability.value == 3.0
    @test fixed_mesh.surface.dwell_time.value == 4.0
    @test fixed_mesh.surface.surface_relaxation.value == 5.0
    @test fixed_mesh.surface.density.value == 6.0

    mesh.size_scale = 0.25
    fixed_mesh_override = fix(mesh)
    @test fixed_mesh_override.geometry isa SizeScaleOverride
    @test size_scale(fixed_mesh_override) == 0.25

    bendy_cylinder = mr.BendyCylinder(
        control_point=[[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        radius=[1.0, 1.0],
        repeats=[2.0, 2.0, 2.0],
        closed=[0, 0, 1],
        spline_order=2,
        nsamples=8,
    )
    fixed_bendy_cylinder = fix(bendy_cylinder)
    @test fixed_bendy_cylinder.geometry isa Repeats.Repeat
    @test first(fixed_bendy_cylinder.geometry.geometry.geometries) isa Meshes.Mesh
end
