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
const Grid = GI.PhysicalGeometries.Grids.Grid
const Transparents = GI.PhysicalGeometries.Transparents
const Transparent = GI.PhysicalGeometries.Transparent
const InsideViews = GI.InsideViews
const IgnoreOverlapping = Transparents.IgnoreOverlapping
const SizeScaleOverride = GI.SizeScaleOverride
const size_scale = GI.size_scale
const has_inside = GI.PhysicalGeometries.has_inside
const has_single_inside = GI.PhysicalGeometries.has_single_inside
const isinside_single = GI.PhysicalGeometries.isinside_single
const inside_indices = GI.PhysicalGeometries.inside_indices
const intersection_type = GI.PhysicalGeometries.intersection_type
const leaf_bounding_boxes = GI.PhysicalGeometries.leaf_bounding_boxes
const has_equivalent_wrapper = GI.PhysicalGeometries.has_equivalent_wrapper
const rewrap = GI.PhysicalGeometries.rewrap
const simplify = GI.PhysicalGeometries.simplify
const bound_intersection_type = GI.PhysicalGeometries.bound_intersection_type
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

@testset "Volume and spin sampling" begin
    fixed_sphere = mr.Geometries.fix(mr.Spheres(radius=1., surface_density=1.))
    bounding_box = mr.BoundingBox(4.)
    positions, inside = GI.volume_sampling(fixed_sphere, bounding_box, 1.)
    @test length(positions) == length(inside)
    @test all(indices isa Vector for indices in inside)

    samples = mr.spin_sampling(fixed_sphere, bounding_box, 1.)
    @test all(typeof(sample).parameters[1] == 0 for sample in samples)
    @test all(sample.isinside isa Vector for sample in samples)

    snapshot = mr.Snapshot(10, bounding_box, fixed_sphere; longitudinal=0.5)
    @test length(snapshot) == 10
    @test all(mr.longitudinal(snapshot) .== 5.0)
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
    @test BoundingBoxes.upper(planar_grid.grid_bounding_box)[3] - BoundingBoxes.lower(planar_grid.grid_bounding_box)[3] == 1.0
    @test size(planar_grid.indices) == (2, 2, 1)

    # Even positive extents smaller than a cell are expanded symmetrically.
    thin_box = BoundingBoxes.InternalBoundingBox([0.125, 1.0, 0.0], [1.0, 2.0, 3.0])
    thin_grid = IntersectionGrid([thin_box]; resolution=1.0)
    @test thin_grid.bounding_box == thin_box
    @test BoundingBoxes.lower(thin_grid.grid_bounding_box) == SVector(0.5, 1.0, 2.5)
    @test BoundingBoxes.upper(thin_grid.grid_bounding_box) == SVector(1.5, 3.0, 3.5)
    @test size(thin_grid.indices) == (1, 2, 1)
    @test thin_grid.inv_resolution == SVector(1.0, 1.0, 1.0)

    # Zero-extent geometry uses a 1 um fallback for automatic and Inf resolution.
    point_box = BoundingBoxes.InternalBoundingBox{3}(0.0, [3.0, 4.0, 5.0])
    for resolution in (nothing, 0.5, Inf)
        point_grid = IntersectionGrid([point_box]; resolution)
        thickness = resolution == 0.5 ? 0.5 : 1.0
        @test point_grid.bounding_box == point_box
        @test BoundingBoxes.lower(point_grid.grid_bounding_box) == SVector(3.0, 4.0, 5.0) .- thickness / 2
        @test BoundingBoxes.upper(point_grid.grid_bounding_box) == SVector(3.0, 4.0, 5.0) .+ thickness / 2
        @test size(point_grid.indices) == (1, 1, 1)
        @test all(isfinite, point_grid.inv_resolution)
        @test point_grid.indices[1, 1, 1] == [1]
    end

    # Inf still means one cell, with finite bounds based on the original extent.
    planar_single = IntersectionGrid([planar_box]; resolution=Inf)
    planar_auto = IntersectionGrid([planar_box])
    for single in (planar_single, planar_auto)
        @test single.bounding_box == planar_box
        @test BoundingBoxes.lower(single.grid_bounding_box) == SVector(-1.0, -1.0, -1.0)
        @test BoundingBoxes.upper(single.grid_bounding_box) == SVector(1.0, 1.0, 1.0)
        @test size(single.indices) == (1, 1, 1)
    end

    # Traversal from outside the geometry stays at a physical, rather than eps, scale.
    intersect = GI.PhysicalGeometries.find_intersection
    wall = Shift(BaseObstructions.InfiniteWall(), [3.0])
    triangle = BaseObstructions.FullTriangle(SVector(0.0, 0.0, 3.0), SVector(1.0, 0.0, 3.0), SVector(0.0, 1.0, 3.0))
    for resolution in (nothing, 0.5, Inf)
        @test intersect(Grid(wall; resolution), SVector(2.0), SVector(4.0)) == intersect(wall, SVector(2.0), SVector(4.0))
        start, destination = SVector(0.1, 0.1, 2.0), SVector(0.1, 0.1, 4.0)
        @test intersect(Grid(triangle; resolution), start, destination) == intersect(triangle, start, destination)
    end
end

@testset "equivalent wrappers" begin
    sphere = BaseObstructions.Sphere(1.0)
    larger = BaseObstructions.Sphere(2.0)
    repeated = Repeat(sphere, [4.0, 4.0, 4.0])
    repeated_larger = Repeat(larger, [4.0, 4.0, 4.0])
    repeated_shifted = Repeat(Shift(sphere, [0.5, 0.0, 0.0]), [4.0, 4.0, 4.0])
    different_period = Repeat(sphere, [4.0, 4.0, 5.0])
    almost_equal_period = Repeat(sphere, [4.0, 4.0, nextfloat(4.0)])

    @test has_equivalent_wrapper(repeated, repeated)
    for other in (repeated_larger, repeated_shifted)
        @test has_equivalent_wrapper(repeated, other)
        @test has_equivalent_wrapper(other, repeated)
    end
    @test has_equivalent_wrapper(repeated_larger, repeated_shifted)
    for other in (different_period, almost_equal_period, sphere,
                  Repeat(BaseObstructions.InfiniteCylinder(1.0), [4.0, 4.0]))
        @test !has_equivalent_wrapper(repeated, other)
        @test !has_equivalent_wrapper(other, repeated)
    end
    # Identical unsupported types still use the conservative fallback.
    for unsupported in (sphere, Shift(sphere, ones(3)), Scale(sphere, 2.0),
                        Rotate(sphere, Matrix{Float64}(I, 3, 3)), IgnoreOverlapping(sphere), Grid(sphere))
        @test !has_equivalent_wrapper(unsupported, unsupported)
    end

    @test has_equivalent_wrapper(GeometryVector([repeated, repeated_larger]))
    @test has_equivalent_wrapper(GeometryTuple((repeated, repeated_larger, repeated_shifted)))
    @test !has_equivalent_wrapper(GeometryVector([repeated, different_period]))
    @test !has_equivalent_wrapper(GeometryTuple((repeated, repeated_shifted, different_period)))
    @test !has_equivalent_wrapper(GeometryTuple((repeated, sphere)))
    @test !has_equivalent_wrapper(GeometryVector([sphere, larger]))
    @test !has_equivalent_wrapper(GeometryVector(BaseObstructions.Sphere[]))
    @test !has_equivalent_wrapper(GeometryTuple{3}(()))
    @test has_equivalent_wrapper(GeometryVector([repeated]))
    @test has_equivalent_wrapper(GeometryTuple((repeated,)))
    @test !has_equivalent_wrapper(GeometryVector([sphere]))
    @test !has_equivalent_wrapper(GeometryTuple((sphere,)))
    group = GeometryTuple((repeated, repeated_larger))
    @test !has_equivalent_wrapper(group, group)
    @test !has_equivalent_wrapper(GeometryTuple((group, group)))
end

@testset "extract common wrapper" begin
    pass = Val(:extract_common_wrapper)
    periods = SVector(6.0, 6.0, 6.0)
    first_child = Shift(BaseObstructions.Sphere(0.5), [-1.2, 0.0, 0.0])
    second_child = Shift(BaseObstructions.Sphere(0.6), [1.3, 0.0, 0.0])
    first_repeat, second_repeat = Repeat(first_child, periods), Repeat(second_child, periods)
    vector = GeometryVector([first_repeat, second_repeat])
    tuple_group = GeometryTuple((first_repeat, second_repeat))
    for group in (vector, tuple_group)
        extracted = simplify(group, pass)
        @test extracted isa Repeat{3}
        @test extracted.repeats == periods
        @test extracted.geometry isa (group isa GeometryVector ? GeometryVector : GeometryTuple)
        @test length(extracted.geometry) == 2
        @test extracted.geometry[1] === first_child
        @test extracted.geometry[2] === second_child
        @test group[1] === first_repeat
        @test group[2] === second_repeat
        @test simplify(extracted.geometry, pass) === extracted.geometry
    end

    # A test-only wrapper lets us exercise consecutive layers without nesting repeats.
    struct ExtractionWrapper{N, P<:PhysicalGeometry{N}} <: Transparent{N, P}
        geometry::P
    end
    GI.PhysicalGeometries.has_equivalent_wrapper(::ExtractionWrapper{N}, ::ExtractionWrapper{N}) where {N} = true
    GI.PhysicalGeometries.rewrap(::ExtractionWrapper{N}, geometry::PhysicalGeometry{N}) where {N} =
        ExtractionWrapper{N, typeof(geometry)}(geometry)
    wrapped_first = ExtractionWrapper{3, typeof(first_repeat)}(first_repeat)
    wrapped_second = ExtractionWrapper{3, typeof(second_repeat)}(second_repeat)
    for group in (GeometryVector([wrapped_first, wrapped_second]), GeometryTuple((wrapped_first, wrapped_second)))
        extracted = simplify(group, pass)
        @test extracted isa ExtractionWrapper
        @test extracted.geometry isa (group isa GeometryVector ? GeometryVector : GeometryTuple)
        @test extracted.geometry[1] === first_repeat
        @test extracted.geometry[2] === second_repeat
        second_pass = simplify(extracted, pass)
        @test second_pass isa ExtractionWrapper
        @test second_pass.geometry isa Repeat
        @test second_pass.geometry.repeats == periods
        inner_group = second_pass.geometry.geometry
        @test inner_group isa (group isa GeometryVector ? GeometryVector : GeometryTuple)
        @test inner_group[1] === first_child
        @test inner_group[2] === second_child
        @test simplify(inner_group, pass) === inner_group
        @test simplify(second_pass, pass) === second_pass
        complete = simplify(group)
        @test typeof(complete) === typeof(second_pass)
        @test complete.geometry.repeats == periods
        @test complete.geometry.geometry.geometries == inner_group.geometries
        @test simplify(complete) === complete
        @test group[1] === wrapped_first
        @test group[2] === wrapped_second
    end

    # Full extraction of each child separately yields different types, whereas
    # one pass allows the parent vector to accept a uniform intermediate form.
    different_repeat = Repeat(second_child, [8., 8., 8.])
    wrapped_different = ExtractionWrapper{3, typeof(different_repeat)}(different_repeat)
    matching = GeometryVector([wrapped_first, wrapped_second])
    differing = GeometryVector([wrapped_first, wrapped_different])
    incremental = GeometryVector([matching, differing])
    greedy_children = map(simplify, incremental.geometries)
    @test typeof(greedy_children[1]) !== typeof(greedy_children[2])
    intermediate = simplify(incremental, pass)
    @test intermediate isa ExtractionWrapper
    @test intermediate.geometry isa GeometryVector
    @test intermediate.geometry[1] isa GeometryVector
    @test intermediate.geometry[2] isa GeometryVector
    @test intermediate.geometry[1].geometries == [first_repeat, second_repeat]
    @test intermediate.geometry[2].geometries == [first_repeat, different_repeat]
    @test simplify(intermediate, pass) === intermediate
    complete = simplify(incremental)
    @test typeof(complete) === typeof(intermediate)
    @test complete.geometry[1].geometries == intermediate.geometry[1].geometries
    @test complete.geometry[2].geometries == intermediate.geometry[2].geometries
    @test simplify(complete) === complete
    @test incremental[1] === matching
    @test incremental[2] === differing

    # Rewrapping uses the new child and rebuilds derived bounds.
    replacement = BaseObstructions.Sphere(1.0)
    replaced = rewrap(first_repeat, replacement)
    @test replaced.geometry === replacement
    @test replaced.repeats == periods
    @test replaced.normalized_bounding_box != first_repeat.normalized_bounding_box
    @test BoundingBoxes.lower(replaced.normalized_bounding_box) == SVector(-1 / 6, -1 / 6, -1 / 6)

    for group in (GeometryVector([first_repeat]), GeometryTuple((first_repeat,)))
        extracted = simplify(group, pass)
        @test extracted isa Repeat
        @test length(extracted.geometry) == 1
        @test extracted.geometry[1] === first_child
    end

    different = Repeat(second_child, [6.0, 6.0, 7.0])
    for unchanged in (
        GeometryVector([first_repeat, different]),
        GeometryTuple((first_repeat, second_repeat, different)),
        GeometryTuple((first_repeat, second_child)),
        GeometryVector([first_child, second_child]),
        GeometryVector(typeof(first_repeat)[]),
        GeometryTuple{3}(()),
        GeometryTuple((Shift(first_repeat, zeros(3)), Shift(second_repeat, zeros(3)))),
    )
        @test simplify(unchanged, pass) === unchanged
    end

    finite = Shift(BaseObstructions.FiniteCylinder([0., 0., -1.], [0., 0., 1.], 0.4), [1.3, 0., 0.])
    heterogeneous = GeometryTuple((first_repeat, Repeat(finite, periods)))
    extracted = simplify(heterogeneous, pass)
    @test extracted.geometry isa GeometryTuple
    @test extracted.geometry.geometries === (first_child, finite)

    # Lifting moves the repeat component before the group component in paths.
    lift_path(path) = (path[2], path[1], path[3:end]...)
    intersect = GI.PhysicalGeometries.find_intersection
    property_index = GI.PhysicalGeometries.to_property_index
    get_params = GI.PhysicalGeometries.get_intersection_params
    properties = Properties.GeometryVectorProperties([10.0, 20.0])
    rng = Random.MersenneTwister(72)
    for original in (vector, tuple_group, heterogeneous)
        lifted = simplify(original, pass)
        for shift in (SVector(0., 0., 0.), SVector(6., -6., 0.))
            start, destination = SVector(-2.5, 0., 0.) + shift, SVector(2.5, 0., 0.) + shift
            hit = intersect(original, start, destination)
            lifted_hit = intersect(lifted, start, destination)
            @test lifted_hit == lift_path(hit)
            @test property_index(original, hit[1:end-2]) == property_index(lifted, lifted_hit[1:end-2])
            @test Properties.get_value(properties, property_index(lifted, lifted_hit[1:end-2])) == 10.0
            @test get_params(original, start, destination, hit).normal ≈ get_params(lifted, start, destination, lifted_hit).normal
            next_hit = intersect(original, start, destination, hit)
            @test intersect(lifted, start, destination, lifted_hit) == lift_path(next_hit)
        end
        for _ in 1:30
            start = SVector{3, Float64}(16 .* rand(rng, 3) .- 8)
            destination = SVector{3, Float64}(16 .* rand(rng, 3) .- 8)
            hit = intersect(original, start, destination)
            @test intersect(lifted, start, destination) == (isnothing(hit) ? nothing : lift_path(hit))
            @test Set(inside_indices(lifted, start)) == Set(lift_path(index) for index in inside_indices(original, start))
        end
        # Probe known interior points, including copies outside the central period.
        for point in (SVector(-1.2, 0., 0.), SVector(1.3, 0., 0.), SVector(7.3, 0., 0.))
            original_inside = inside_indices(original, point)
            lifted_inside = inside_indices(lifted, point)
            @test !isempty(original_inside)
            @test Set(lifted_inside) == Set(lift_path(index) for index in original_inside)
            @test Set(property_index(original, index) for index in original_inside) ==
                Set(property_index(lifted, index) for index in lifted_inside)
        end
    end

    # Surface suffixes survive extraction and gap flags are preserved.
    gaps = Shift(BaseObstructions.FiniteCylinder([0., 0., -1.], [0., 0., 1.], 0.4; caps_are_gaps=true), [1.3, 0., 0.])
    caps = GeometryTuple((first_repeat, Repeat(gaps, periods)))
    lifted_caps = simplify(caps, pass)
    start, destination = SVector(1.3, 0., -2.), SVector(1.3, 0., 2.)
    cap_hit = intersect(caps, start, destination)
    lifted_cap_hit = intersect(lifted_caps, start, destination)
    @test lifted_cap_hit == lift_path(cap_hit)
    @test cap_hit[3] == 2
    @test get_params(lifted_caps, start, destination, lifted_cap_hit).hit_gap
    @test property_index(caps, cap_hit[1:end-2]) == property_index(lifted_caps, lifted_cap_hit[1:end-2])

    walls = GeometryVector([
        Repeat(Shift(BaseObstructions.InfiniteWall(), [-1.]), [6.]),
        Repeat(Shift(BaseObstructions.InfiniteWall(), [1.]), [6.]),
    ])
    @test simplify(walls, pass) isa Repeat{1}
    @test intersect(simplify(walls, pass), SVector(-2.), SVector(2.)) ==
        lift_path(intersect(walls, SVector(-2.), SVector(2.)))
end

@testset "simplify geometry" begin
    pass = Val(:extract_common_wrapper)
    periods = SVector(6., 6., 6.)
    a = Shift(BaseObstructions.Sphere(0.3), [-1., 0., 0.])
    b = Shift(BaseObstructions.Sphere(0.4), [1., 0., 0.])
    c = Shift(BaseObstructions.Sphere(0.3), [0., 1., 0.])
    d = Shift(BaseObstructions.Sphere(0.4), [0., -1., 0.])
    first = GeometryVector([Repeat(a, periods), Repeat(b, periods)])
    second = GeometryVector([Repeat(c, periods), Repeat(d, periods)])
    mixed = GeometryVector([Repeat(a, periods), Repeat(b, [8., 8., 8.])])

    # Children simplify to different types: reject every proposed vector change.
    inconsistent = GeometryVector([first, mixed])
    @test simplify(first, pass) isa Repeat
    @test simplify(mixed, pass) === mixed
    @test simplify(inconsistent, pass) === inconsistent
    @test inconsistent[1] === first
    @test inconsistent[2] === mixed
    @test first[1].geometry === a

    # A uniform element-type change is accepted without converting to a tuple.
    other_periods = GeometryVector([Repeat(c, [8., 8., 8.]), Repeat(d, [8., 8., 8.])])
    uniform = GeometryVector([first, other_periods])
    uniform_result = simplify(uniform, pass)
    @test uniform_result isa GeometryVector
    @test uniform_result !== uniform
    @test eltype(uniform_result) !== eltype(uniform)
    @test isconcretetype(eltype(uniform_result))
    @test all(child -> child isa Repeat, uniform_result)
    @test uniform_result[1].geometry.geometries == [a, b]
    @test uniform_result[2].geometry.geometries == [c, d]
    @test simplify(uniform_result, pass) === uniform_result

    # Tuple groups accept heterogeneous child changes while preserving their family.
    tuple_input = GeometryTuple((first, mixed))
    tuple_result = simplify(tuple_input, pass)
    @test tuple_result isa GeometryTuple
    @test tuple_result[1] isa Repeat
    @test tuple_result[2] === mixed
    @test tuple_input[1] === first
    @test simplify(tuple_result, pass) === tuple_result

    # Child extraction exposes equivalent repeats for the parent to extract.
    for original in (GeometryVector([first, second]), GeometryTuple((first, second)))
        result = simplify(original, pass)
        @test result isa Repeat
        @test result.repeats == periods
        @test result.geometry isa (original isa GeometryVector ? GeometryVector : GeometryTuple)
        @test result.geometry[1].geometries == [a, b]
        @test result.geometry[2].geometries == [c, d]
        @test original[1] === first
        @test original[2] === second
        @test simplify(result, pass) === result

        # The repeat index moves ahead of both original group indices.
        move_repeat(path) = (path[3], path[1], path[2], path[4:end]...)
        intersect = GI.PhysicalGeometries.find_intersection
        property_index = GI.PhysicalGeometries.to_property_index
        for shift in (SVector(0., 0., 0.), SVector(6., -6., 0.))
            start, destination = SVector(-2., 0., 0.) + shift, SVector(2., 0., 0.) + shift
            hit = intersect(original, start, destination)
            new_hit = intersect(result, start, destination)
            @test new_hit == move_repeat(hit)
            @test intersect(result, start, destination, new_hit) ==
                move_repeat(intersect(original, start, destination, hit))
            @test property_index(original, hit[1:end-2]) == property_index(result, new_hit[1:end-2])
            @test GI.PhysicalGeometries.get_intersection_params(original, start, destination, hit).normal ≈
                GI.PhysicalGeometries.get_intersection_params(result, start, destination, new_hit).normal
            for point in (a.shift + shift, b.shift + shift, c.shift + shift, d.shift + shift)
                before = inside_indices(original, point)
                after = inside_indices(result, point)
                @test !isempty(before)
                @test Set(after) == Set(move_repeat(index) for index in before)
                @test Set(property_index(original, index) for index in before) ==
                    Set(property_index(result, index) for index in after)
            end
        end
    end

    rotation = [0. -1. 0.; 1. 0. 0.; 0. 0. 1.]
    wrapped = Shift(Scale(Rotate(IgnoreOverlapping(SizeScaleOverride(first, 0.1)), rotation), 2.), [1., 2., 3.])
    wrapped_result = simplify(wrapped, pass)
    @test wrapped_result isa Shift
    @test wrapped_result.shift == wrapped.shift
    @test wrapped_result.geometry isa Scale
    @test wrapped_result.geometry.scale == 2.
    @test wrapped_result.geometry.geometry isa Rotate
    @test wrapped_result.geometry.geometry.matrix == rotation
    @test wrapped_result.geometry.geometry.geometry isa IgnoreOverlapping
    override = wrapped_result.geometry.geometry.geometry.geometry
    @test override isa SizeScaleOverride
    @test override.size_scale == 0.1
    @test override.geometry isa Repeat
    @test wrapped.geometry.geometry.geometry.geometry.geometry === first
    @test simplify(wrapped_result, pass) === wrapped_result

    cylinders = GeometryVector([
        Repeat(Shift(BaseObstructions.InfiniteCylinder(0.3), [-1., 0.]), [6., 6.]),
        Repeat(Shift(BaseObstructions.InfiniteCylinder(0.3), [1., 0.]), [6., 6.]),
    ])
    projected = Rotate(cylinders, Matrix{Float64}(I, 3, 2))
    projected_result = simplify(projected, pass)
    @test projected_result isa Rotate{3, 2}
    @test projected_result.geometry isa Repeat{2}
    @test projected_result.matrix == projected.matrix

    # Repeat also simplifies its child, rebuilding derived state when it changes.
    struct SimplifiableSphere <: PhysicalGeometry{3}
        radius::Float64
    end
    BoundingBoxes.InternalBoundingBox(geometry::SimplifiableSphere) = BoundingBoxes.InternalBoundingBox{3}(geometry.radius)
    GI.PhysicalGeometries.simplify(geometry::SimplifiableSphere, ::Val) = geometry
    GI.PhysicalGeometries.simplify(geometry::SimplifiableSphere, ::Val{:extract_common_wrapper}) =
        geometry.radius < 0.25 ? geometry : BaseObstructions.Sphere(geometry.radius)
    repeated = Repeat(SimplifiableSphere(0.5), periods)
    repeated_result = simplify(repeated, pass)
    @test repeated_result isa Repeat{3, BaseObstructions.Sphere}
    @test repeated_result.repeats == repeated.repeats
    @test repeated_result.normalized_bounding_box == repeated.normalized_bounding_box
    @test repeated.geometry isa SimplifiableSphere

    # Default traversal forwards other passes without applying the extraction rule.
    other_pass = Val(:traverse_only)
    for original in (first, inconsistent, uniform, tuple_input, wrapped, projected, repeated)
        @test simplify(original, other_pass) === original
    end
    struct PassProbe <: PhysicalGeometry{3}
        calls::Vector{Symbol}
    end
    function GI.PhysicalGeometries.simplify(probe::PassProbe, ::Val{S}) where {S}
        push!(probe.calls, S)
        probe
    end
    probes = [PassProbe(Symbol[]), PassProbe(Symbol[])]
    probe_tree = Shift(IgnoreOverlapping(GeometryTuple((GeometryVector(probes), probes[1]))), zeros(3))
    @test simplify(probe_tree, other_pass) === probe_tree
    @test probes[1].calls == [:traverse_only, :traverse_only]
    @test probes[2].calls == [:traverse_only]

    # A single extraction pass visits each original leaf only once.
    foreach(probe -> empty!(probe.calls), probes)
    wrapped_probes = GeometryVector([ExtractionWrapper{3, PassProbe}(probe) for probe in probes])
    @test simplify(wrapped_probes, pass) isa ExtractionWrapper
    @test probes[1].calls == [:extract_common_wrapper]
    @test probes[2].calls == [:extract_common_wrapper]

    # A different local rule can change descendants without extracting groups.
    GI.PhysicalGeometries.simplify(geometry::SimplifiableSphere, ::Val{:replace_spheres}) = BaseObstructions.Sphere(geometry.radius)
    probe_repeats = GeometryVector([repeated, Repeat(SimplifiableSphere(0.7), periods)])
    replacements = simplify(probe_repeats, Val(:replace_spheres))
    @test replacements isa GeometryVector
    @test all(child -> child isa Repeat{3, BaseObstructions.Sphere}, replacements)
    @test all(child -> child.geometry isa SimplifiableSphere, probe_repeats)
    @test simplify(replacements, Val(:replace_spheres)) === replacements

    # Reject child changes, but still extract the original equivalent wrappers.
    rejected = GeometryVector([Repeat(SimplifiableSphere(0.1), periods), repeated])
    @test has_equivalent_wrapper(rejected)
    rejected_result = simplify(rejected, pass)
    @test rejected_result isa Repeat
    @test rejected_result.repeats == periods
    @test rejected_result.geometry isa GeometryVector
    @test eltype(rejected_result.geometry) === SimplifiableSphere
    @test rejected_result.geometry[1] === rejected[1].geometry
    @test rejected_result.geometry[2] === rejected[2].geometry
    @test simplify(rejected_result, pass) === rejected_result
    @test simplify(rejected, other_pass) === rejected
    @test simplify(rejected[1], pass) === rejected[1]
    @test simplify(rejected[2], pass) isa Repeat{3, BaseObstructions.Sphere}

    sphere = BaseObstructions.Sphere(0.5)
    @test simplify(sphere) === sphere
    for unchanged in (sphere, first[1], GeometryVector([a, b]), GeometryTuple((a, b)),
                      GeometryVector(typeof(a)[]), GeometryTuple{3}(()), Shift(sphere, ones(3)),
                      Grid(GeometryVector([a, b]); resolution=0.5))
        @test simplify(unchanged, pass) === unchanged
        @test simplify(unchanged) === unchanged
    end

    struct UnsupportedSimplification <: PhysicalGeometry{3} end
    unsupported = UnsupportedSimplification()
    @test_throws MethodError simplify(unsupported)
    @test_throws MethodError simplify(unsupported, pass)
    @test_throws MethodError simplify(unsupported, other_pass)
    @test_throws MethodError simplify(GeometryVector([unsupported]), pass)
    @test_throws MethodError simplify(GeometryTuple((sphere, unsupported)), pass)
    @test_throws MethodError simplify(Shift(unsupported, zeros(3)), pass)
    for primitive in (BaseObstructions.InfiniteWall(), BaseObstructions.InfiniteCylinder(1.),
                      BaseObstructions.FiniteCylinder([0., 0., 0.], [0., 0., 1.], 0.5))
        @test simplify(primitive, pass) === primitive
    end
end

@testset "simplify liminal geometry" begin
    Liminal = GI.PhysicalGeometries.LiminalGeometries
    sphere = BaseObstructions.Sphere(1.)
    smaller = BaseObstructions.Sphere(0.75)
    singleton = GeometryTuple((sphere,))
    geometry = Liminal.FixedLiminalGeometry(
        PhysicalGeometry{3}[singleton, smaller], [1., 3.], 0.2;
        tortuosity_strength=0.7,
    )
    @test simplify(geometry, Val(:traverse_only)) === geometry
    @test simplify(geometry) === geometry

    # Remove a group layer so the old cached surface paths are no longer valid.
    GI.PhysicalGeometries.simplify(group::GeometryTuple{3, Tuple{BaseObstructions.Sphere}}, ::Val{:unwrap_liminal_test}) =
        group[1].radius >= 1. ? group[1] : group
    sampling = geometry.outer_surface_sampling
    @test length(sampling.positions) > 1
    sampling.index_to_sample = 2
    original_positions = copy(sampling.positions)
    original_indices = copy(sampling.surface_indices)
    result = simplify(geometry, Val(:unwrap_liminal_test))
    @test result isa Liminal.FixedLiminalGeometry
    @test result !== geometry
    @test result.geometries == [sphere, smaller]
    @test geometry.geometries[1] === singleton
    @test result.number_fractions == geometry.number_fractions
    @test result.extracellular_fraction == geometry.extracellular_fraction
    @test result.tortuosity_strength == geometry.tortuosity_strength
    @test result.total_surface_area == geometry.total_surface_area
    @test result.weighted_cell_volume == geometry.weighted_cell_volume
    @test result.spherical_surface_area === geometry.spherical_surface_area
    @test result.normal_second_moment_tensor == geometry.normal_second_moment_tensor
    @test result.tortuosity_transform == geometry.tortuosity_transform
    @test result.inverse_tortuosity_transform == geometry.inverse_tortuosity_transform
    fresh = result.outer_surface_sampling
    @test fresh !== sampling
    @test fresh.lock !== sampling.lock
    @test fresh.index_to_sample == 1
    @test !isempty(fresh.positions)
    @test eltype(fresh.surface_indices) == Tuple{Bool}
    @test all(length(index) == 1 for index in fresh.surface_indices)
    @test all(zip(fresh.positions, fresh.surface_indices, fresh.cell_indices)) do (position, index, cell)
        params = GI.PhysicalGeometries.get_intersection_params(
            result.geometries[cell], position, position, (index..., 0.),
        )
        norm(params.normal) ≈ 1.
    end
    @test sampling.positions == original_positions
    @test sampling.surface_indices == original_indices
    @test sampling.index_to_sample == 2
    @test simplify(result, Val(:unwrap_liminal_test)) === result
    @test simplify(result) === result

    # Liminal collections can accept heterogeneous child types after rewriting.
    homogeneous = rewrap(geometry, [singleton, GeometryTuple((smaller,))])
    heterogeneous = simplify(homogeneous, Val(:unwrap_liminal_test))
    @test heterogeneous.geometries[1] === sphere
    @test heterogeneous.geometries[2] === homogeneous.geometries[2]
    @test length(Base.uniontypes(eltype(heterogeneous.geometries))) == 2
    @test heterogeneous.weighted_cell_volume == geometry.weighted_cell_volume
    @test heterogeneous.spherical_surface_area === geometry.spherical_surface_area

    # Partial precomputation skips volume estimation; complete precomputation
    # also skips area and large spherical-statistics sampling during rewrap.
    struct LiminalEstimationProbe <: Transparent{3, BaseObstructions.Sphere}
        geometry::BaseObstructions.Sphere
        volume_calls::Base.RefValue{Int}
        surface_calls::Base.RefValue{Int}
        sampling_calls::Base.RefValue{Int}
    end
    function GI.PhysicalGeometries.LiminalGeometries._estimate_inside_volume(probe::LiminalEstimationProbe; kwargs...)
        probe.volume_calls[] += 1
        4π / 3 * probe.geometry.radius^3
    end
    function GI.PhysicalGeometries.estimate_surface(probe::LiminalEstimationProbe; kwargs...)
        probe.surface_calls[] += 1
        GI.PhysicalGeometries.estimate_surface(probe.geometry; kwargs...)
    end
    function GI.PhysicalGeometries.random_surface_positions(
        probe::LiminalEstimationProbe, density::Properties.GeometryProperties,
        box::BoundingBoxes.InternalBoundingBox{3}, scale_density; kwargs...,
    )
        probe.sampling_calls[] += 1
        GI.PhysicalGeometries.random_surface_positions(probe.geometry, density, box, scale_density; kwargs...)
    end
    probe = LiminalEstimationProbe(sphere, Ref(0), Ref(0), Ref(0))
    volume = 4π / 3
    precomputed_volume = Liminal.FixedLiminalGeometry(
        [probe], [1.], 0.2; precomputed=(weighted_cell_volume=volume,),
    )
    @test precomputed_volume.weighted_cell_volume == volume
    @test probe.volume_calls[] == 0
    @test probe.surface_calls[] == 1
    @test probe.sampling_calls[] == 2
    rebuilt = rewrap(precomputed_volume, [probe])
    @test probe.volume_calls[] == 0
    @test probe.surface_calls[] == 1
    @test probe.sampling_calls[] == 3
    @test rebuilt.spherical_surface_area === precomputed_volume.spherical_surface_area
    @test rebuilt.outer_surface_sampling !== precomputed_volume.outer_surface_sampling
    @test rebuilt.outer_surface_sampling.index_to_sample == 1
end

@testset "leaf bounding boxes" begin
    sphere = BaseObstructions.Sphere(1.0)
    cylinder = BaseObstructions.InfiniteCylinder(2.0)
    wall = BaseObstructions.InfiniteWall()
    triangle = BaseObstructions.FullTriangle(SVector(0., 0., 0.), SVector(1., 0., 0.), SVector(0., 1., 0.))
    for primitive in (sphere, cylinder, wall, triangle)
        @test leaf_bounding_boxes(primitive) == BoundingBoxes.InternalBoundingBox(primitive)
    end

    # The tuple order matches the side/first-cap/second-cap surface identifiers.
    finite = BaseObstructions.FiniteCylinder([1., 2., 3.], [1., 2., 7.], 1., 2.)
    finite_boxes = leaf_bounding_boxes(finite)
    @test finite_boxes isa Tuple{Vararg{BoundingBoxes.InternalBoundingBox{3}, 3}}
    @test finite_boxes[1] == BoundingBoxes.InternalBoundingBox(finite)
    @test BoundingBoxes.lower(finite_boxes[2]) == SVector(0., 1., 3.)
    @test BoundingBoxes.upper(finite_boxes[2]) == SVector(2., 3., 3.)
    @test BoundingBoxes.lower(finite_boxes[3]) == SVector(-1., 0., 7.)
    @test BoundingBoxes.upper(finite_boxes[3]) == SVector(3., 4., 7.)
    gaps = BaseObstructions.FiniteCylinder([1., 2., 3.], [1., 2., 7.], 1., 2.; caps_are_gaps=true)
    @test leaf_bounding_boxes(gaps) == finite_boxes
    tilted = BaseObstructions.FiniteCylinder([0., 0., 0.], [3., 4., 0.], 1., 2.)
    tilted_boxes = leaf_bounding_boxes(tilted)
    @test BoundingBoxes.lower(tilted_boxes[2]) ≈ SVector(-0.8, -0.6, -1.)
    @test BoundingBoxes.upper(tilted_boxes[2]) ≈ SVector(0.8, 0.6, 1.)
    @test BoundingBoxes.lower(tilted_boxes[3]) ≈ SVector(1.4, 2.8, -2.)
    @test BoundingBoxes.upper(tilted_boxes[3]) ≈ SVector(4.6, 5.2, 2.)

    nested = GeometryTuple((sphere, GeometryVector([finite])))
    boxes = leaf_bounding_boxes(nested)
    @test boxes isa Tuple
    @test boxes[1] == BoundingBoxes.InternalBoundingBox(sphere)
    @test boxes[2] isa Vector
    @test boxes[2] == [finite_boxes]
    for group in (GeometryVector([sphere]), GeometryVector([sphere]; bounding_box=true), GeometryVector([sphere]; grid=true))
        @test leaf_bounding_boxes(group) == [BoundingBoxes.InternalBoundingBox(sphere)]
    end
    @test leaf_bounding_boxes(GeometryTuple{3}(())) === ()
    @test isempty(leaf_bounding_boxes(GeometryVector(BaseObstructions.Sphere[])))

    # Transform every box while retaining nested vectors and component tuples.
    rotation = [0. -1. 0.; 1. 0. 0.; 0. 0. 1.]
    transformed = Shift(Rotate(Scale(nested, 2.), rotation), [10., 20., 30.])
    transformed_boxes = leaf_bounding_boxes(transformed)
    @test transformed_boxes isa Tuple
    @test transformed_boxes[2] isa Vector
    @test transformed_boxes[2][1] isa Tuple
    @test transformed_boxes[1] == BoundingBoxes.InternalBoundingBox{3}(2., [10., 20., 30.])
    @test BoundingBoxes.lower(transformed_boxes[2][1][2]) == SVector(4., 20., 36.)
    @test BoundingBoxes.upper(transformed_boxes[2][1][2]) == SVector(8., 24., 36.)
    @test BoundingBoxes.lower(transformed_boxes[2][1][3]) == SVector(2., 18., 44.)
    @test BoundingBoxes.upper(transformed_boxes[2][1][3]) == SVector(10., 26., 44.)
    @test leaf_bounding_boxes(SizeScaleOverride(IgnoreOverlapping(nested), 0.5)) == boxes
    @test isempty(leaf_bounding_boxes(Shift(GeometryVector(BaseObstructions.Sphere[]), ones(3))))

    # Triangle boxes include the mesh's automatically generated gap closure.
    mesh = Mesh(
        [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1, 2, 3], [1, 2, 4], [1, 3, 4]],
    )
    mesh_boxes = leaf_bounding_boxes(mesh)
    @test mesh_boxes isa Vector{BoundingBoxes.InternalBoundingBox{3}}
    @test length(mesh_boxes) == length(mesh.indices)
    @test mesh.first_index_of_gap <= length(mesh_boxes)
    for (index, vertices) in enumerate(mesh.indices)
        points = mesh.vertices[vertices]
        @test BoundingBoxes.lower(mesh_boxes[index]) == reduce((a, b) -> min.(a, b), points)
        @test BoundingBoxes.upper(mesh_boxes[index]) == reduce((a, b) -> max.(a, b), points)
    end

    repeat = Repeat(sphere, [4., 4., 4.])
    @test_throws PublicBoundingBoxes.BoundingBoxNotSupported leaf_bounding_boxes(repeat)
    @test_throws PublicBoundingBoxes.BoundingBoxNotSupported leaf_bounding_boxes(GeometryTuple((sphere, repeat)))
    projected = Rotate(cylinder, Matrix{Float64}(I, 3, 2))
    @test_throws PublicBoundingBoxes.BoundingBoxNotSupported leaf_bounding_boxes(projected)
end

@testset "single intersection candidate" begin
    intersect = GI.PhysicalGeometries.find_intersection
    sphere = BaseObstructions.Sphere(1.)
    start, destination = SVector(-5., 0., 0.), SVector(5., 0., 0.)
    children = [Shift(sphere, [-2., 0., 0.]), Shift(sphere, [2., 0., 0.])]
    for group in (GeometryVector(children), GeometryVector(children; bounding_box=true),
                  GeometryVector(children; grid=true, grid_resolution=1.), GeometryTuple(Tuple(children)))
        @test intersect(group, start, destination, nothing, nothing, nothing) == intersect(group, start, destination)
        @test intersect(group, start, destination, nothing, nothing, (1,)) == (1, false, 0.2)
        @test intersect(group, start, destination, nothing, nothing, (2,)) == (2, false, 0.6)
        @test intersect(group, start, destination, (1, false, 0.2), nothing, (2,)) == (2, false, 0.6)
        @test intersect(group, start, destination, (2, false, 0.6), nothing, (2,)) === nothing
        @test intersect(group, start .+ SVector(0., 4., 0.), destination .+ SVector(0., 4., 0.), nothing, nothing, (2,)) === nothing
        @test_throws ArgumentError intersect(group, start, destination, nothing, nothing, ())
        @test_throws BoundsError intersect(group, start, destination, nothing, nothing, (3,))
        @test_throws ArgumentError intersect(group, start, destination, nothing, nothing, (1, 1))
    end

    nested = GeometryTuple((sphere, GeometryVector(children)))
    wrapped = SizeScaleOverride(IgnoreOverlapping(nested), 0.25)
    @test intersect(wrapped, start, destination, nothing, nothing, (2, 2)) == (2, 2, false, 0.6)
    rotation = [0. -1. 0.; 1. 0. 0.; 0. 0. 1.]
    transformed = Shift(Rotate(Scale(wrapped, 2.), rotation), [1., 2., 3.])
    a, b = SVector(1., -8., 3.), SVector(1., 12., 3.)
    @test intersect(transformed, a, b, nothing, nothing, (2, 2)) == (2, 2, false, 0.6)

    for (primitive, a, b) in (
        (sphere, start, destination),
        (BaseObstructions.InfiniteCylinder(1.), SVector(-2., 0.), SVector(2., 0.)),
        (BaseObstructions.InfiniteWall(), SVector(-1.), SVector(1.)),
        (BaseObstructions.FullTriangle(SVector(0., 0., 0.), SVector(1., 0., 0.), SVector(0., 1., 0.)),
         SVector(0.1, 0.1, -1.), SVector(0.1, 0.1, 1.)),
    )
        hit = intersect(primitive, a, b)
        @test hit !== nothing
        @test intersect(primitive, a, b, nothing, nothing, ()) == hit
        @test intersect(primitive, a, b, hit, nothing, ()) == intersect(primitive, a, b, hit)
        @test_throws ArgumentError intersect(primitive, a, b, nothing, nothing, (1,))
    end

    finite = BaseObstructions.FiniteCylinder([0., 0., -1.], [0., 0., 1.], 1.)
    a, b = SVector(0., 0., -2.), SVector(0., 0., 2.)
    @test intersect(finite, a, b) == (2, false, 0.25)
    @test intersect(finite, a, b, nothing, nothing, (1,)) === nothing
    @test intersect(finite, a, b, nothing, nothing, (2,)) == (2, false, 0.25)
    # A selected later cap must not be discarded in favor of the overall nearest hit.
    @test intersect(finite, a, b, nothing, nothing, (3,)) == (3, false, 0.75)
    side_a, side_b = SVector(-2., 0., 0.), SVector(2., 0., 0.)
    @test intersect(finite, side_a, side_b, nothing, nothing, (1,)) == (1, false, 0.25)
    @test intersect(finite, side_a, side_b, nothing, nothing, (2,)) === nothing
    @test intersect(finite, side_a, side_b, nothing, nothing, (3,)) === nothing
    previous = (2, true, 0.25)
    @test intersect(finite, a, b, previous, nothing, (2,)) === nothing
    @test intersect(finite, a, b, previous, nothing, (3,)) == (3, true, 0.75)
    @test intersect(finite, a, b, (2, false, 0.25), nothing, (3,)) === nothing
    @test intersect(GeometryTuple((sphere, finite)), a, b, nothing, nothing, (2, 3)) == (2, 3, false, 0.75)
    gaps = BaseObstructions.FiniteCylinder([0., 0., -1.], [0., 0., 1.], 1.; caps_are_gaps=true)
    cap_hit = intersect(gaps, a, b, nothing, nothing, (3,))
    @test GI.PhysicalGeometries.get_intersection_params(gaps, a, b, cap_hit).hit_gap
    for path in ((), (0,), (4,), (1, 2))
        @test_throws ArgumentError intersect(finite, a, b, nothing, nothing, path)
    end

    mesh = Mesh(
        [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1, 2, 3], [1, 2, 4], [1, 3, 4]],
    )
    a, b = SVector(-1., 0.1, 0.1), SVector(2., 0.1, 0.1)
    for index in eachindex(mesh.indices)
        triangle_hit = intersect(Meshes.triangle(mesh, index), a, b)
        expected = isnothing(triangle_hit) ? nothing : (index, triangle_hit...)
        @test intersect(mesh, a, b, nothing, nothing, (index,)) == expected
        @test intersect(mesh, a, b, (index, false, 0.), nothing, (index,)) === nothing
    end
    gap_index = mesh.first_index_of_gap
    gap_hit = intersect(mesh, a, b, nothing, nothing, (gap_index,))
    @test gap_hit !== nothing
    @test gap_hit[end] > intersect(mesh, a, b)[end]
    @test GI.PhysicalGeometries.get_intersection_params(mesh, a, b, gap_hit).hit_gap
    @test intersect(mesh, a, b, (0, false, 0.), nothing, (gap_index,)) == gap_hit
    @test intersect(GeometryTuple((sphere, mesh)), a, b, nothing, nothing, (2, gap_index)) == (2, gap_hit...)
    @test_throws ArgumentError intersect(mesh, a, b, nothing, nothing, ())
    @test_throws ArgumentError intersect(mesh, a, b, nothing, nothing, (1, 1))
    @test_throws BoundsError intersect(mesh, a, b, nothing, nothing, (length(mesh.indices) + 1,))

    # Instrumented children verify direct dispatch and child-local state routing.
    struct CandidateProbe <: PhysicalGeometry{3}
        calls::Vector{Any}
    end
    GI.PhysicalGeometries.find_intersection_requires_inside(::Type{CandidateProbe}) = Val(true)
    function GI.PhysicalGeometries.find_intersection(
        probe::CandidateProbe, a::SVector{3, Float64}, b::SVector{3, Float64},
        previous=nothing, inside=nothing, candidate=nothing,
    )
        push!(probe.calls, (; previous, inside=isnothing(inside) ? nothing : collect(inside), candidate))
        (false, 0.75)
    end
    skipped, first_probe, selected = CandidateProbe(Any[]), CandidateProbe(Any[]), CandidateProbe(Any[])
    probes = IgnoreOverlapping(GeometryTuple((skipped, GeometryVector([first_probe, selected]))))
    cached = InsideViews.InsideView([(1,), (2, 1), (2, 2)])
    previous = (2, 2, true, 0.25)
    @test intersect(probes, start, destination, previous, cached, (2, 2)) == (2, 2, false, 0.75)
    @test isempty(skipped.calls)
    @test isempty(first_probe.calls)
    @test only(selected.calls) == (; previous=(true, 0.25), inside=[()], candidate=())
    empty!(selected.calls)
    @test intersect(probes, start, destination, (1, false, 0.2), cached, (2, 2)) == (2, 2, false, 0.75)
    @test only(selected.calls).previous === nothing

    repeat = Repeat(sphere, [4., 4., 4.])
    @test intersect(repeat, start, destination, nothing, nothing, nothing) == intersect(repeat, start, destination)
    @test_throws ArgumentError intersect(repeat, start, destination, nothing, nothing, ())
    @test_throws ArgumentError intersect(repeat, start, destination, nothing, nothing, (1,))
end

@testset "transparent intersection grid" begin
    intersect = GI.PhysicalGeometries.find_intersection
    sphere = BaseObstructions.Sphere(1.)
    finite = BaseObstructions.FiniteCylinder([0., 0., -1.], [0., 0., 1.], 0.5)
    tree = GeometryTuple((GeometryVector([
        Shift(sphere, [-2., 0., 0.]), Shift(sphere, [2., 0., 0.]),
    ]), finite))
    grid = Grid(tree; resolution=0.5)
    @test grid isa Transparent{3}
    @test grid.geometry === tree
    @test grid.indices == [(1, 1), (1, 2), (2, 1), (2, 2), (2, 3)]
    @test eltype(grid.indices) == intersection_type(typeof(tree))
    @test BoundingBoxes.InternalBoundingBox(grid) == BoundingBoxes.InternalBoundingBox(tree)
    @test leaf_bounding_boxes(grid) == leaf_bounding_boxes(tree)
    @test size_scale(grid) == size_scale(tree)
    @test Grid(sphere).indices == [()]

    rng = Random.MersenneTwister(19)
    for resolution in (nothing, 0.5, Inf)
        accelerated = Grid(tree; resolution)
        for _ in 1:40
            start = SVector{3, Float64}(10 .* rand(rng, 3) .- 5)
            destination = SVector{3, Float64}(10 .* rand(rng, 3) .- 5)
            @test intersect(accelerated, start, destination) == intersect(tree, start, destination)
            @test inside_indices(accelerated, start) == inside_indices(tree, start)
        end
    end

    start, destination = SVector(-5., 0., 0.), SVector(5., 0., 0.)
    hit = intersect(grid, start, destination)
    @test hit == intersect(tree, start, destination)
    @test intersect(grid, start, destination, hit) == intersect(tree, start, destination, hit)
    @test intersect(grid, start, start) === nothing
    @test intersect(grid, SVector(-5., 4., 0.), SVector(5., 4., 0.)) === nothing
    @test GI.PhysicalGeometries.to_property_index(grid, hit[1:(end - 2)]) == (1, 1)
    @test GI.PhysicalGeometries.to_inside_index(grid, (2, 3)) == (2,)
    params = GI.PhysicalGeometries.get_intersection_params(grid, start, destination, hit)
    @test params.normal ≈ SVector(-1., 0., 0.)
    @test !params.hit_gap
    for (a, b) in ((SVector(0., 0., -2.), SVector(0., 0., 2.)),
                   (SVector(-0.75, 0., 0.), SVector(0.75, 0., 0.)))
        @test intersect(grid, a, b) == intersect(tree, a, b)
    end

    # Cached inside states must retain other obstructions for overlap handling.
    overlapping = IgnoreOverlapping(GeometryTuple((
        sphere, Shift(sphere, [0.5, 0., 0.]),
    )))
    overlap_grid = Grid(overlapping; resolution=0.5)
    a, b = SVector(0.25, 0., 0.), SVector(2., 0., 0.)
    cached = InsideViews.InsideView([(1,), (2,)])
    overlap_hit = intersect(overlap_grid, a, b, nothing, cached)
    @test overlap_hit == intersect(overlapping, a, b, nothing, cached)
    @test GI.PhysicalGeometries.get_intersection_params(overlap_grid, a, b, overlap_hit, cached).hit_gap

    # A broad box can yield a farther hit than a candidate visited later.
    broad = Shift(BaseObstructions.Sphere(5.), [0., 4., 0.])
    near = Shift(BaseObstructions.Sphere(0.25), [-4., 0., 0.])
    competing = GeometryTuple((broad, near))
    @test intersect(Grid(competing; resolution=0.5), start, destination) == intersect(competing, start, destination)
    @test intersect(Grid(competing; resolution=0.5), start, destination)[1] == 2
    # Equal-distance hits use tree order even if the second box is visited first.
    tied = GeometryTuple((Shift(sphere, [-2., 0., 0.]), broad))
    @test intersect(Grid(tied; resolution=0.5), start, destination) == intersect(tied, start, destination)

    mesh = Mesh(
        [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1, 2, 3], [1, 2, 4], [1, 3, 4]],
    )
    mesh_grid = Grid(mesh; resolution=0.25)
    @test mesh_grid.indices == [(index,) for index in eachindex(mesh.indices)]
    a, b = SVector(0.1, 0.1, 0.1), SVector(1., 1., 1.)
    mesh_hit = intersect(mesh_grid, a, b)
    @test mesh_hit == intersect(mesh, a, b)
    @test mesh_hit[1] >= mesh.first_index_of_gap
    @test GI.PhysicalGeometries.get_intersection_params(mesh_grid, a, b, mesh_hit).hit_gap
    @test isinside_single(mesh_grid, a) == isinside_single(mesh, a)
    @test intersect(mesh_grid, a, b, mesh_hit) == intersect(mesh, a, b, mesh_hit)

    rotation = [0. -1. 0.; 1. 0. 0.; 0. 0. 1.]
    transformed = Shift(Rotate(Scale(tree, 2.), rotation), [1., 2., 3.])
    @test intersect(Grid(transformed; resolution=1.), SVector(1., -8., 3.), SVector(1., 12., 3.)) ==
        intersect(transformed, SVector(1., -8., 3.), SVector(1., 12., 3.))
    cylinder = BaseObstructions.InfiniteCylinder(1.)
    @test intersect(Grid(cylinder), SVector(-2., 0.), SVector(2., 0.)) == intersect(cylinder, SVector(-2., 0.), SVector(2., 0.))
    walls = GeometryVector([Shift(BaseObstructions.InfiniteWall(), [-1.]), Shift(BaseObstructions.InfiniteWall(), [1.])])
    @test intersect(Grid(walls; resolution=0.5), SVector(-2.), SVector(2.)) == intersect(walls, SVector(-2.), SVector(2.))

    periodic = Repeat(tree, [8., 8., 8.])
    periodic_grid = Repeat(grid, [8., 8., 8.])
    @test periodic_grid.geometry === grid
    @test intersect(periodic_grid, start, destination) == intersect(periodic, start, destination)
    for shift in (SVector(8., 0., 0.), SVector(-8., 8., 0.))
        a, b = start + shift, destination + shift
        previous = intersect(periodic, a, b)
        @test intersect(periodic_grid, a, b) == previous
        @test intersect(periodic_grid, a, b, previous) == intersect(periodic, a, b, previous)
        @test inside_indices(periodic_grid, shift) == inside_indices(periodic, shift)
    end
    @test_throws ArgumentError intersect(grid, start, destination, nothing, nothing, ())
    @test_throws ArgumentError Grid(GeometryTuple{3}(()))
    @test_throws ArgumentError Grid(sphere; resolution=0.)
    @test_throws PublicBoundingBoxes.BoundingBoxNotSupported Grid(periodic)

    # A leaf spanning several cells is still checked only once by GridIterator.
    struct GridProbe <: BaseObstructions.BaseObstruction{3}
        calls::Base.RefValue{Int}
    end
    BoundingBoxes.InternalBoundingBox(::GridProbe) = BoundingBoxes.InternalBoundingBox{3}(1.)
    function GI.PhysicalGeometries.find_intersection(
        probe::GridProbe, a::SVector{3, Float64}, b::SVector{3, Float64},
        previous=nothing, inside=nothing, candidate=nothing,
    )
        @test candidate == ()
        probe.calls[] += 1
        nothing
    end
    selected, skipped = GridProbe(Ref(0)), GridProbe(Ref(0))
    probes = GeometryTuple((selected, Shift(skipped, [0., 4., 0.])))
    @test intersect(Grid(probes; resolution=0.25), start, destination) === nothing
    @test selected.calls[] == 1
    @test skipped.calls[] == 0
end

@testset "intersection types" begin
    sphere_type = BaseObstructions.Sphere

    @test intersection_type(sphere_type) == Tuple{}
    @test intersection_type(Mesh) == Tuple{Int}
    @test intersection_type(GeometryVector{3, sphere_type}) == Tuple{Int}
    @test intersection_type(Repeat{3, sphere_type}) == Tuple{SVector{3, Int}}
    @test intersection_type(Shift{3, sphere_type}) == Tuple{}
    @test intersection_type(SizeScaleOverride{3, sphere_type}) == Tuple{}
    @test intersection_type(GeometryTuple{3, Tuple{BaseObstructions.Sphere, Mesh}}) == Union{Tuple{Int}, Tuple{Int, Int}}
end

@testset "bound intersection types" begin
    sphere = BaseObstructions.Sphere(1.)
    zero_density = Properties.GeometryLeafProperties(0.)
    density = Properties.GeometryLeafProperties(1.)

    @test bound_intersection_type(sphere, zero_density) == Union{}
    @test bound_intersection_type(sphere, density) == Tuple{}

    sphere_group = GeometryVector([sphere, sphere])
    group_density = Properties.GeometryVectorProperties([0., 1.])
    @test bound_intersection_type(sphere_group, group_density) == Tuple{Int}
    @test bound_intersection_type(sphere_group, Properties.GeometryVectorProperties([0., 0.])) == Union{}

    repeated_sphere = Repeat(sphere, [4., 4., 4.])
    @test bound_intersection_type(repeated_sphere, density) == Tuple{SVector{3, Int}}
    @test bound_intersection_type(Shift(sphere, zero(SVector{3})), density) == Tuple{}

    nested = GeometryTuple((sphere, GeometryTuple((sphere, sphere))))
    nested_density_a = Properties.GeometryTupleProperties((1., Properties.GeometryTupleProperties((0., 0.))))
    nested_density_b = Properties.GeometryTupleProperties((0., Properties.GeometryTupleProperties((1., 0.))))
    nested_group = GeometryVector([nested, nested])
    nested_densities = Properties.GeometryVectorProperties([nested_density_a, nested_density_b])
    @test bound_intersection_type(nested_group, nested_densities) == Union{Tuple{Int, Int}, Tuple{Int, Int, Int}}

    fixed_geometry = mr.fix(mr.Spheres(radius=1., surface_density=1., dwell_time=1.))
    @test GI.bound_intersection_type(fixed_geometry) ==
        GI.Intersection{Tuple{Int}, Tuple{Int}}
    @test mr.Reflections.possible_reflection_types(fixed_geometry) ==
        Union{Nothing, mr.Reflections.Reflection{GI.Intersection{Tuple{Int}, Tuple{Int}}}}
    zero_fixed_geometry = mr.fix(mr.Spheres(radius=1.))
    @test mr.Reflections.possible_reflection_types(zero_fixed_geometry) === Nothing
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

    free = nothing
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
    @test !has_inside(typeof(BaseObstructions.InfiniteWall()))
    @test !has_inside(typeof(BaseObstructions.FullTriangle(
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
    )))

    @testset "Inside views" begin
        indices = [(1,), (2, 1), (2, 2)]
        view = InsideViews.InsideView(indices)
        @test collect(view) == indices

        first_child = InsideViews.child_view(view, 1)
        second_child = InsideViews.child_view(view, 2)
        @test collect(first_child) == [()]
        @test collect(second_child) == [(1,), (2,)]
        @test !InsideViews.has_other(first_child, ())
        @test InsideViews.has_other(second_child, (1,))

        properties = Properties.GeometryTupleProperties((1., Properties.GeometryTupleProperties((2., 3.))))
        @test Properties.get_value(properties, view) == 6.
        @test Properties.get_value(properties.properties[2], second_child) == 5.

        liminal_view = InsideViews.InsideView([
            ((1, :offset_a), 1),
            ((2, :offset_b), 1),
        ])
        @test collect(InsideViews.child_view(liminal_view, (1, :offset_a))) == [(1,)]

        narrowed_view = InsideViews.InsideView(indices, 2, 3, 0)
        deeper_view = InsideViews.InsideView(
            narrowed_view.indices,
            narrowed_view.first,
            narrowed_view.last,
            narrowed_view.depth + 1,
        )
        @test (deeper_view.first, deeper_view.last, deeper_view.depth) == (2, 3, 1)
        @test collect(deeper_view) == [(1,), (2,)]
    end

    @testset "Spherical surface area" begin
        Random.seed!(1234)
        N = 10_000
        normals = [
            let normal = SVector(randn(), randn(), randn())
                normal / norm(normal)
            end
            for _ in 1:N
        ]
        sample_weight = 4π / N
        spherical_area = GI.PhysicalGeometries.LiminalGeometries.SphericalSurfaceArea(
            normals;
            degree=8,
            sample_weight=sample_weight,
        )
        @test length(spherical_area.coefficients) == 81
        for direction in (SVector(1., 0., 0.), SVector(1., 2., 3.))
            unit_direction = direction / norm(direction)
            expected = sample_weight * sum(abs(normal ⋅ unit_direction) for normal in normals) / 2
            @test GI.PhysicalGeometries.LiminalGeometries.projected_surface_area(
                spherical_area,
                direction,
            ) ≈ expected atol=0.01
        end
        closed_normals = reduce(vcat, [[normal, -normal] for normal in normals])
        closed_surface_area = GI.PhysicalGeometries.LiminalGeometries.SphericalSurfaceArea(
            closed_normals;
            degree=8,
            sample_weight=sample_weight / 2,
        )
        for direction in (SVector(1., 0., 0.), SVector(1., 2., 3.))
            unit_direction = direction / norm(direction)
            weights = max.(0., -([normal ⋅ unit_direction for normal in closed_normals]))
            expected = sum(
                weight * (normal ⋅ unit_direction)^2
                for (weight, normal) in zip(weights, closed_normals)
            ) / (sum(weights))
            @test GI.PhysicalGeometries.LiminalGeometries.normal_second_moment(
                closed_surface_area,
                direction,
            ) ≈ expected atol=0.01
        end
        @test_throws ArgumentError GI.PhysicalGeometries.LiminalGeometries.projected_surface_area(
            spherical_area,
            zero(SVector{3}),
        )
        @test_throws ArgumentError GI.PhysicalGeometries.LiminalGeometries.normal_second_moment(
            closed_surface_area,
            zero(SVector{3}),
        )

        triangle_normals = [
            SA[1., 0., 0.],
            SA[-cos(π / 3), sin(2π / 3), 0.],
            SA[-cos(π / 3), -sin(2π / 3), 0.],
        ]
        triangle_surface_area = GI.PhysicalGeometries.LiminalGeometries.SphericalSurfaceArea(
            triangle_normals;
            degree=32,
            sample_weight=1.,
        )
        for direction in (
            SA[1., 0., 0.],
            SA[cos(π / 12), sin(π / 12), 0.],
            SA[cos(π / 6), sin(π / 6), 0.],
            SA[-1., 0., 0.],
            SA[0., 0., 1.],
            SA[cos(π / 6), 0., sin(π / 6)],
        )
            cosines = [normal ⋅ direction for normal in triangle_normals]
            incoming_weights = max.(0., -cosines)
            expected_projected_area = sum(abs, cosines) / 2
            expected_moment = if sum(incoming_weights) == 0
                0.
            else
                sum(
                    weight * cosine^2
                    for (weight, cosine) in zip(incoming_weights, cosines)
                ) / sum(incoming_weights)
            end

            @test GI.PhysicalGeometries.LiminalGeometries.projected_surface_area(
                triangle_surface_area,
                direction,
            ) ≈ expected_projected_area atol=0.03
            @test GI.PhysicalGeometries.LiminalGeometries.normal_second_moment(
                triangle_surface_area,
                direction,
            ) ≈ expected_moment atol=0.03
        end
    end

    @testset "Intersection inside requirements" begin
        sphere = BaseObstructions.Sphere(1.0)
        requires_inside = GI.PhysicalGeometries.get_intersection_params_requires_inside
        @test (@inferred requires_inside(typeof(sphere))) === Val(false)
        @test (@inferred requires_inside(typeof(IgnoreOverlapping(sphere)))) === Val(true)
        @test (@inferred requires_inside(typeof(Shift(IgnoreOverlapping(sphere), zero(SVector{3}))))) === Val(true)
        @test (@inferred requires_inside(typeof(GeometryTuple{3}(())))) === Val(false)
        @test (@inferred requires_inside(typeof(GeometryTuple((sphere, sphere))))) === Val(false)
        @test (@inferred requires_inside(typeof(GeometryTuple((sphere, IgnoreOverlapping(sphere)))))) === Val(true)
        nested = GeometryTuple((sphere, Shift(GeometryTuple((sphere, IgnoreOverlapping(sphere))), zero(SVector{3}))))
        @test (@inferred requires_inside(typeof(nested))) === Val(true)
    end

    @testset "Inside index conversion" begin
        sphere = BaseObstructions.Sphere(1.0)
        triangle = BaseObstructions.FullTriangle(
            SVector(0., 0., 0.), SVector(1., 0., 0.), SVector(0., 1., 0.),
        )
        mesh = Mesh(
            [SVector(0., 0., 0.), SVector(1., 0., 0.), SVector(0., 1., 0.), SVector(0., 0., 1.)],
            [SVector(1, 2, 3), SVector(1, 2, 4), SVector(1, 3, 4)],
        )
        to_inside = GI.PhysicalGeometries.to_inside_index
        @test to_inside(sphere, ()) == ()
        @test to_inside(triangle, ()) === nothing
        @test to_inside(mesh, (2,)) == ()
        @test to_inside(GeometryVector([mesh]), (1, 2)) == (1,)
        repeated_mesh = Repeat(mesh, [2., 2., 2.])
        @test to_inside(repeated_mesh, (SVector(0, 0, 0), 2)) == (SVector(0, 0, 0),)
    end

    @testset "Finite cylinders" begin
        cylinder = BaseObstructions.FiniteCylinder(
            SVector(0., 0., 0.), SVector(0., 0., 2.), 1.,
        )
        variable = BaseObstructions.FiniteCylinder(
            SVector(0., 0., 0.), SVector(0., 0., 2.), 1., 2.,
        )
        keyword_variable = BaseObstructions.FiniteCylinder(
            SVector(0., 0., 0.), SVector(0., 0., 2.), 1.; radius_second=2.,
        )
        @test keyword_variable == variable
        @test_throws ArgumentError BaseObstructions.FiniteCylinder(
            SVector(0., 0., 0.), SVector(0., 0., 0.), 1.,
        )
        @test_throws ArgumentError BaseObstructions.FiniteCylinder(
            SVector(0., 0., 0.), SVector(0., 0., 2.), 0.,
        )
        @test has_inside(typeof(cylinder))
        @test has_single_inside(typeof(cylinder))
        @test intersection_type(typeof(cylinder)) == Tuple{Int}
        @test isinside_single(cylinder, SVector(0., 0., 1.))
        @test !isinside_single(cylinder, SVector(1.1, 0., 1.))
        @test isinside_single(cylinder, SVector(1.1, 0., 1.), (true, 0.5))
        @test BoundingBoxes.lower(BoundingBoxes.InternalBoundingBox(cylinder)) == SVector(-1., -1., 0.)
        @test BoundingBoxes.upper(BoundingBoxes.InternalBoundingBox(cylinder)) == SVector(1., 1., 2.)

        side_hit = GI.PhysicalGeometries.find_intersection(
            cylinder, SVector(2., 0., 1.), SVector(-2., 0., 1.),
        )
        @test side_hit == (1, false, 0.25)
        inside_side_hit = GI.PhysicalGeometries.find_intersection(
            cylinder, SVector(0., 0., 1.), SVector(2., 0., 1.),
        )
        @test inside_side_hit == (1, true, 0.5)
        cap_hit = GI.PhysicalGeometries.find_intersection(
            cylinder, SVector(0., 0., 1.), SVector(0., 0., 3.),
        )
        @test cap_hit == (3, true, 0.5)
        @test isnothing(GI.PhysicalGeometries.find_intersection(
            cylinder, SVector(0., 0., 2.), SVector(0., 0., 3.), (3, true, 0.5),
        ))
        variable_hit = GI.PhysicalGeometries.find_intersection(
            variable, SVector(2., 0., 1.), SVector(-3., 0., 1.),
        )
        @test variable_hit == (1, false, 0.1)
        @test isnothing(GI.PhysicalGeometries.find_intersection(
            cylinder, SVector(2., 0., 3.), SVector(-2., 0., 3.),
        ))

        side_params = GI.PhysicalGeometries.get_intersection_params(
            cylinder, SVector(2., 0., 1.), SVector(-2., 0., 1.), side_hit, nothing,
        )
        @test side_params.normal == SVector(1., 0., 0.)
        gap_cylinder = BaseObstructions.FiniteCylinder(
            SVector(0., 0., 0.), SVector(0., 0., 2.), 1.; caps_are_gaps=true,
        )
        gap_params = GI.PhysicalGeometries.get_intersection_params(
            gap_cylinder, SVector(0., 0., 1.), SVector(0., 0., 3.), cap_hit, nothing,
        )
        @test gap_params.normal == SVector(0., 0., -1.)
        @test gap_params.hit_gap
        @test !GI.PhysicalGeometries.get_intersection_params(
            cylinder, SVector(0., 0., 1.), SVector(0., 0., 3.), cap_hit, nothing,
        ).hit_gap
        gap_surface = GI.PhysicalGeometries.estimate_surface(
            gap_cylinder; density=100., minimum_samples=1_000,
        )
        gap_surface_included = GI.PhysicalGeometries.estimate_surface(
            gap_cylinder; density=100., minimum_samples=1_000, include_gap=true,
        )
        @test gap_surface.area ≈ 4π rtol=0.1
        @test gap_surface_included.area ≈ 6π rtol=0.1
        gap_density = Properties.GeometryLeafProperties(1.)
        Random.seed!(1234)
        _, samples_without_gap = GI.PhysicalGeometries.random_surface_positions(
            gap_cylinder, gap_density, BoundingBoxes.InternalBoundingBox(gap_cylinder), 100.;
            include_gap=false,
        )
        Random.seed!(1234)
        _, samples_with_gap = GI.PhysicalGeometries.random_surface_positions(
            gap_cylinder, gap_density, BoundingBoxes.InternalBoundingBox(gap_cylinder), 100.;
            include_gap=true,
        )
        @test all(index[1] == 1 for index in samples_without_gap)
        @test any(index[1] != 1 for index in samples_with_gap)

        overlapping_cylinder = IgnoreOverlapping(cylinder)
        overlapping_params = GI.PhysicalGeometries.get_intersection_params(
            overlapping_cylinder,
            SVector(0., 0., 1.),
            SVector(0., 0., 3.),
            cap_hit,
            [()],
        )
        @test !overlapping_params.hit_gap

        cylinder_geometry = mr.FiniteCylinders(
            position=[[0., 0., 0.], [0., 0., 2.]],
            radius=[1., 1.],
            connected_to=[0, 1],
            use_spherical_endpoint=[false, false],
        )
        template = mr.fix(cylinder_geometry)
        fixed_cylinder = GI.FixedGeometry(
            cylinder, template.volume, template.surface, template.susceptibility,
        )
        side_detection = GI.detect_intersection(
            fixed_cylinder, SVector(0., 0., 1.), SVector(2., 0., 1.),
        )
        @test side_detection.distance ≈ 0.5
        @test side_detection.inside
        cap_detection = GI.detect_intersection(
            fixed_cylinder, SVector(0., 0., 1.), SVector(0., 0., 3.),
        )
        @test cap_detection.distance ≈ 0.5
        @test cap_detection.inside
        @test !cap_detection.hit_gap

        @test length(GI.geometry_mesh(variable; nsamples=8)[1].vertices) == 18
        @test length(GI.geometry_mesh(variable; nsamples=8)[1].triangles) == 32
    end

    @testset "Monte Carlo geometry measures" begin
        sphere = BaseObstructions.Sphere(1.)
        @testset "inside_sampling" begin
            sample_density = 600.
            approximately_uniform(values) = begin
                bins = 10
                all(
                    abs(count(value -> (bin - 1) / bins <= value < bin / bins, values) /
                        length(values) - 1 / bins) < 0.025
                    for bin in 1:bins
                )
            end

            Random.seed!(1234)
            sphere_positions, sphere_indices = GI.inside_sampling(
                sphere, sample_density;
                rng=Random.MersenneTwister(1234),
            )
            @test length(sphere_positions) / sample_density ≈ 4π / 3 rtol=0.04
            @test all(index == [()] for index in sphere_indices)
            @test mean(sum(abs2, position) for position in sphere_positions) ≈ 3 / 5 rtol=0.02
            @test approximately_uniform([norm(position)^3 for position in sphere_positions])
            @test approximately_uniform([
                (position[3] / norm(position) + 1) / 2
                for position in sphere_positions
            ])
            @test approximately_uniform([
                mod(atan(position[2], position[1]), 2π) / (2π)
                for position in sphere_positions
            ])

            circle = BaseObstructions.InfiniteCylinder(1.)
            circle_positions, circle_indices = GI.PhysicalGeometries.inside_sampling(
                circle, sample_density; rng=Random.MersenneTwister(1234),
            )
            @test length(circle_positions) / sample_density ≈ π rtol=0.04
            @test all(index == [()] for index in circle_indices)
            @test approximately_uniform([sum(abs2, position) for position in circle_positions])
            @test approximately_uniform([
                mod(atan(position[2], position[1]), 2π) / (2π)
                for position in circle_positions
            ])

            tapered_cylinder = BaseObstructions.FiniteCylinder(
                SVector(0., 0., 0.), SVector(0., 0., 3.), 1., 2.,
            )
            cylinder_positions, cylinder_indices = GI.PhysicalGeometries.inside_sampling(
                tapered_cylinder, sample_density;
                rng=Random.MersenneTwister(2345),
            )
            @test length(cylinder_positions) / sample_density ≈ 7π rtol=0.04
            @test all(index == [()] for index in cylinder_indices)
            @test mean(position[3] for position in cylinder_positions) ≈ 3 * 17 / 28 rtol=0.015
            axial_cdf = [
                let fraction = position[3] / tapered_cylinder.length
                    (fraction + fraction^2 + fraction^3 / 3) / (7 / 3)
                end
                for position in cylinder_positions
            ]
            normalized_radius_squared = [
                (position[1]^2 + position[2]^2) /
                    (1 + position[3] / tapered_cylinder.length)^2
                for position in cylinder_positions
            ]
            radial_angles = [
                mod(atan(-position[1], position[2]), 2π) / (2π)
                for position in cylinder_positions
            ]
            @test approximately_uniform(axial_cdf)
            @test approximately_uniform(normalized_radius_squared)
            @test approximately_uniform(radial_angles)

            overlapping_spheres = mr.fix(mr.Spheres(
                position=[[0., 0., 0.], [1., 0., 0.]],
                radius=1.,
                overlapping=true,
            ))
            Random.seed!(3456)
            overlap_positions, overlap_indices = GI.inside_sampling(
                overlapping_spheres, 1500.,
            )
            union_volume = 9π / 4
            @test length(overlap_positions) / 1500 ≈ union_volume rtol=0.035
            @test all(!isempty(index) for index in overlap_indices)
            overlap_fraction = count(index -> length(index) == 2, overlap_indices) /
                length(overlap_indices)
            @test overlap_fraction ≈ (5π / 12) / union_volume rtol=0.06

            mesh = Mesh(
                [SVector(0., 0., 0.), SVector(1., 0., 0.), SVector(0., 1., 0.), SVector(0., 0., 1.)],
                [SVector(1, 2, 3), SVector(1, 2, 4), SVector(1, 3, 4)],
            )
            mesh_fixed = GI.FixedGeometry(
                mesh,
                Properties.GeometryLeafProperties(0.),
                Properties.GeometryLeafProperties(0.),
                nothing,
            )
            Random.seed!(4567)
            mesh_positions, mesh_indices = GI.inside_sampling(mesh_fixed, 1000.)
            @test length(mesh_positions) / 1000 ≈ 1 / 6 rtol=0.1
            @test all(index == [()] for index in mesh_indices)
        end
        volume = GI.PhysicalGeometries.estimate_volume(
            sphere; nsamples=100_000, rng=Random.MersenneTwister(1234),
        )
        @test volume ≈ 4π / 3 atol=0.2
        @test GI.PhysicalGeometries.estimate_volume(
            sphere; nsamples=10, minimum_samples=1, maximum_attempts=2,
            rng=Random.MersenneTwister(1234),
        ) isa Float64
        @test_throws ArgumentError GI.PhysicalGeometries.estimate_volume(sphere; minimum_samples=-1)
        @test_throws ArgumentError GI.PhysicalGeometries.estimate_volume(sphere; maximum_attempts=0)

        Random.seed!(1234)
        surface = GI.PhysicalGeometries.estimate_surface(
            sphere; density=100., minimum_samples=1_000,
        )
        @test surface.area ≈ 4π rtol=0.1
        @test length(surface.positions) == length(surface.normals)
        @test length(surface.positions) == length(surface.indices)
        @test length(surface.positions) == length(surface.weights)

        cylinder = BaseObstructions.FiniteCylinder(
            SVector(0., 0., 0.), SVector(0., 0., 2.), 1.,
        )
        volume = GI.PhysicalGeometries.estimate_volume(
            cylinder; nsamples=100_000, rng=Random.MersenneTwister(1234),
        )
        @test volume ≈ 2π atol=0.2

        Random.seed!(1234)
        surface = GI.PhysicalGeometries.estimate_surface(
            cylinder; density=100., minimum_samples=1_000,
        )
        @test surface.area ≈ 6π rtol=0.1

        overlapping = GeometryTuple((
            Shift(BaseObstructions.Sphere(1.), SVector(0., 0., 0.)),
            Shift(BaseObstructions.Sphere(1.), SVector(0.5, 0., 0.)),
        ))
        Random.seed!(1234)
        all_surface = GI.PhysicalGeometries.estimate_surface(
            overlapping; density=100., minimum_samples=1_000,
        )
        Random.seed!(1234)
        outer_surface = GI.PhysicalGeometries.estimate_surface(
            overlapping; density=100., minimum_samples=1_000, outer=true,
        )
        @test outer_surface.area < all_surface.area
        @test length(outer_surface.positions) < length(all_surface.positions)

        fixed_sphere = mr.fix(mr.Spheres(radius=1.))
        @test mr.volume(fixed_sphere; nsamples=100_000, rng=Random.MersenneTwister(1234)) ≈ 4π / 3 atol=0.2
        Random.seed!(1234)
        @test mr.surface(fixed_sphere; density=100., minimum_samples=1_000) ≈ 4π rtol=0.1

        annulus = mr.fix(mr.Annuli(inner=0.5, outer=1.))
        annulus_box = mr.BoundingBox(2.)
        Random.seed!(1234)
        all_annulus_surface = mr.surface(
            annulus; bounding_box=annulus_box, density=100., minimum_samples=1_000,
        )
        Random.seed!(1234)
        outer_annulus_surface = mr.surface(
            annulus; bounding_box=annulus_box, density=100., minimum_samples=1_000, outer=true,
        )
        axial_length = 2 * annulus_box.upper[3]
        @test all_annulus_surface ≈ axial_length * 2π * (0.5 + 1.) rtol=0.1
        @test outer_annulus_surface ≈ axial_length * 2π rtol=0.1

        repeating_sphere = mr.fix(mr.Spheres(radius=1., repeats=[3., 3., 3.]))
        @test_throws ArgumentError mr.volume(repeating_sphere)
        @test_throws ArgumentError mr.surface(repeating_sphere)

        cells = mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1.)), (2., mr.fix(mr.Spheres(radius=2.)))],
            extracellular_fraction=0.2,
        )
        @test cells.geometries[1][2] isa mr.Spheres
        @test cells.geometries[2][2] isa GI.FixedGeometry
        fixed_cells = mr.fix(cells)
        @test fixed_cells isa GI.FixedGeometry
        @test fixed_cells.geometry isa GI.PhysicalGeometries.LiminalGeometries.FixedLiminalGeometry
        @test_throws ArgumentError GI.PhysicalGeometries.find_intersection(
            fixed_cells.geometry, SVector(0., 0., 0.), SVector(1., 0., 0.), nothing, nothing, (1,),
        )
        sparse_centers = [[1.5 * index, 0., 0.] for index in -10:10]
        sparse_cells = mr.fix(mr.Spheres(
            position=sparse_centers,
            radius=0.1,
            overlapping=true,
        ))
        sparse_liminal = mr.fix(mr.LiminalGeometry(
            geometries=[(1., sparse_cells)],
            extracellular_fraction=0.2,
        ))
        @test sparse_liminal.geometry.weighted_cell_volume ≈
            length(sparse_centers) * 4π * 0.1^3 / 3 rtol=0.04
        geometry_mesh_error = try
            GI.geometry_mesh(fixed_cells.geometry)
            nothing
        catch exception
            exception
        end
        @test geometry_mesh_error isa ArgumentError
        @test geometry_mesh_error.msg == "geometry meshes are not supported for liminal geometries"
        distance_error = try
            GI.distance_to_surface(fixed_cells.geometry, SVector(0., 0., 0.))
            nothing
        catch exception
            exception
        end
        @test distance_error isa ArgumentError
        @test distance_error.msg == "distance_to_surface is not defined for liminal geometries"
        @test GI.size_scale(fixed_cells.geometry) == 1.
        liminal_box = mr.BoundingBox(2.)
        @test GI.PhysicalGeometries.estimate_volume(
            fixed_cells.geometry; bounding_box=liminal_box,
        ) == 0.8 * 4^3
        @test GI.PhysicalGeometries.estimate_surface(
            fixed_cells.geometry; bounding_box=liminal_box,
        ).area ≈ 27/17 * 4^3 * 0.8 rtol=0.1
        @test_throws mr.BoundingBoxNotSupported GI.PhysicalGeometries.estimate_volume(
            fixed_cells.geometry,
        )
        @test_throws mr.BoundingBoxNotSupported GI.PhysicalGeometries.estimate_surface(
            fixed_cells.geometry,
        )
        @test fixed_cells.geometry.outer_surface_sampling isa
            GI.PhysicalGeometries.LiminalGeometries.OuterSurfaceSampling
        @test fixed_cells.geometry.total_surface_area > 0
        spherical_surface_area = fixed_cells.geometry.spherical_surface_area
        spherical_coefficients = copy(spherical_surface_area.coefficients)
        sampling = fixed_cells.geometry.outer_surface_sampling
        GI.PhysicalGeometries.LiminalGeometries.sample!(sampling, fixed_cells.geometry, 100)
        @test !isempty(sampling.positions)
        @test length(sampling.positions) == length(sampling.normals) ==
            length(sampling.cell_indices) == length(sampling.surface_indices)
        @test sampling.weight > 0
        direction = SVector(1., 2., 3.)
        @test GI.projected_surface_area(sampling, direction) ≈
            GI.projected_surface_area(sampling, -direction)
        @test GI.inverse_mean_free_path(fixed_cells.geometry, direction) ≈
            mr.inverse_mean_free_path(fixed_cells, direction)
        @test spherical_surface_area.coefficients == spherical_coefficients
        @test GI.inverse_mean_free_path(fixed_cells.geometry, direction) ≈
            (1 - fixed_cells.geometry.extracellular_fraction) /
            (fixed_cells.geometry.extracellular_fraction * fixed_cells.geometry.weighted_cell_volume) *
            GI.projected_surface_area(spherical_surface_area, direction)
        @test_throws ArgumentError GI.projected_surface_area(sampling, zero(SVector{3}))
        @test fixed_cells.geometry.geometries isa Vector{<:GI.PhysicalGeometry}
        @test eltype(fixed_cells.geometry.geometries) !== GI.PhysicalGeometry
        @test fixed_cells.geometry.number_fractions ≈ [1 / 3, 2 / 3]
        @test sum(fixed_cells.geometry.number_fractions) ≈ 1.
        @test fixed_cells.geometry.extracellular_fraction == 0.2

        sphere = BaseObstructions.Sphere(1.)
        cylinder = Rotate(
            BaseObstructions.InfiniteCylinder(1.),
            [1. 0.; 0. 1.; 0. 0.],
        )
        mixed_cells = GI.PhysicalGeometries.LiminalGeometries.FixedLiminalGeometry(
            [sphere, cylinder], [0.5, 0.5], 0.2,
        )
        @test mixed_cells.total_surface_area ≈ 3π atol=0.1
        @test mixed_cells.weighted_cell_volume ≈ 7π / 6 atol=0.05

        parallel_direction = SVector(0., 0., 1.)
        perpendicular_direction = SVector(1., 0., 0.)
        @test GI.inverse_mean_free_path(mixed_cells, parallel_direction) <
            GI.inverse_mean_free_path(mixed_cells, perpendicular_direction)

        Random.seed!(1234)
        parallel_cells = Int[]
        for _ in 1:100
            hit = GI.PhysicalGeometries.find_intersection(
                mixed_cells,
                SVector(0., 0., -10.),
                SVector(0., 0., 10.),
            )
            hit !== nothing && push!(parallel_cells, hit[1][1])
        end
        @test !isempty(parallel_cells)
        @test all(==(1), parallel_cells)

        Random.seed!(1234)
        perpendicular_cells = Int[]
        for _ in 1:100
            hit = GI.PhysicalGeometries.find_intersection(
                mixed_cells,
                SVector(-10., 0., 0.),
                SVector(10., 0., 0.),
            )
            hit !== nothing && push!(perpendicular_cells, hit[1][1])
        end
        @test 1 in perpendicular_cells
        @test 2 in perpendicular_cells

        sphere_sampling = GI.PhysicalGeometries.LiminalGeometries.OuterSurfaceSampling(
            SVector{3, Float64}[],
            SVector{3, Float64}[],
            Int[],
            Tuple[],
            0.,
            1,
            ReentrantLock(),
        )
        sphere_geometry = GI.PhysicalGeometries.LiminalGeometries.FixedLiminalGeometry(
            [BaseObstructions.Sphere(1.)], [1.], 0.2,
        )
        GI.PhysicalGeometries.LiminalGeometries.sample!(sphere_sampling, sphere_geometry, 100)
        @test all(
            position ⋅ normal > 0
            for (position, normal) in zip(sphere_sampling.positions, sphere_sampling.normals)
        )
        Random.seed!(1234)
        liminal_intersection = GI.PhysicalGeometries.find_intersection(
            sphere_geometry,
            SVector(0., 0., 0.),
            SVector(10., 0., 0.),
        )
        @test liminal_intersection !== nothing
        @test liminal_intersection[end - 1] === false
        @test 0 < liminal_intersection[end] < 1
        @test GI.PhysicalGeometries.find_intersection(
            sphere_geometry,
            SVector(0., 0., 0.),
            SVector(10., 0., 0.),
            liminal_intersection,
        ) !== nothing
        @test fixed_cells.volume.R1 isa GI.Properties.GeometryVectorProperties
        @test length(fixed_cells.volume.R1.properties) == 2
        @test GI.PhysicalGeometries.child_type(typeof(fixed_cells.geometry)) ===
            Core.apply_type(Transformations.Shift, 3, eltype(fixed_cells.geometry.geometries))
        @test GI.PhysicalGeometries.has_inside(typeof(fixed_cells.geometry))
        @test !GI.PhysicalGeometries.has_single_inside(typeof(fixed_cells.geometry))
        error = try
            GI.PhysicalGeometries.inside_indices(fixed_cells.geometry, SVector(0., 0., 0.))
            nothing
        catch exception
            exception
        end
        @test error isa ArgumentError
        @test occursin("use the cached inside indices instead", error.msg)
        @test GI.PhysicalGeometries.inside_indices_eltype(typeof(fixed_cells.geometry)) ==
            Tuple{Tuple{Int, SVector{3, Float64}}, Int}
        @test GI.PhysicalGeometries.intersection_type(typeof(fixed_cells.geometry)) ==
            Tuple{Tuple{Int, SVector{3, Float64}}, Int}
        sampling_box = mr.BoundingBox(5.)
        volume_positions, volume_indices = GI.volume_sampling(fixed_cells, sampling_box, 0.2)
        @test length(volume_positions) == length(volume_indices)
        @test count(isempty, volume_indices) ≈
            0.2 * 0.2 * 10^3 rtol=0.25
        @test count(!isempty, volume_indices) ≈
            0.8 * 0.2 * 10^3 rtol=0.25
        @test all(all(position .>= -5) && all(position .<= 5) for position in volume_positions)
        @test all(
            isempty(index) || all(length(cell_index) == 2 for cell_index in index)
            for index in volume_indices
        )

        surface_cells = mr.fix(mr.LiminalGeometry(
            geometries=[
                (1., mr.Spheres(radius=1., surface_density=1.)),
                (2., mr.Spheres(radius=2., surface_density=1.)),
            ],
            extracellular_fraction=0.2,
        ))
        surface_positions, surface_intersections = GI.random_surface_positions(
            surface_cells, sampling_box, 1.,
        )
        @test length(surface_positions) == length(surface_intersections)
        @test all(all(position .>= -5) && all(position .<= 5) for position in surface_positions)
        @test all(length(intersection.indices) == 2 for intersection in surface_intersections)
        @test all(length(intersection.indices[1][2]) == 3 for intersection in surface_intersections)
        surface_samples = mr.spin_sampling(surface_cells, sampling_box, 0.2)
        @test all(sample isa mr.Spin for sample in surface_samples)
        @test all(
            sample.isinside == mr.isinside(surface_cells, sample.position, sample.reflection.intersection).inside_of
            for sample in surface_samples if !isnothing(sample.reflection)
        )

        Random.seed!(1234)
        annulus = mr.Geometries.fix(mr.Annuli(
            inner=0.5,
            outer=1.,
            inner_surface_density=1.,
            outer_surface_density=1.,
        ))
        annulus_samples = filter(
            sample -> !isnothing(sample.reflection),
            mr.spin_sampling(annulus, mr.BoundingBox(4.), 1.),
        )
        inner_samples = filter(
            sample -> sample.reflection.intersection.indices[1] == 1,
            annulus_samples,
        )
        @test !isempty(inner_samples)
        @test all(
            sample.isinside == mr.isinside(
                annulus,
                sample.position,
                sample.reflection.intersection,
            ).inside_of
            for sample in inner_samples
        )
        @test all(any(index -> index[1] == 2, sample.isinside) for sample in inner_samples)
        @test all(
            any(index -> index == (1, 1), sample.isinside) == sample.reflection.intersection.inside
            for sample in inner_samples
        )
        @test_throws UndefKeywordError mr.LiminalGeometry(geometries=[])
        @test_throws ArgumentError mr.LiminalGeometry(
            geometries=[], extracellular_fraction=0.,
        )
        @test_throws ArgumentError mr.LiminalGeometry(geometries=[(0., mr.Spheres(radius=1.))], extracellular_fraction=0.)
        @test_throws ArgumentError mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1.))], extracellular_fraction=1.1,
        )
        @test_throws ArgumentError mr.fix(mr.LiminalGeometry(
            geometries=[(1., repeating_sphere)], extracellular_fraction=0.,
        ))
        default_tortuosity = mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1.))],
            extracellular_fraction=0.2,
        )
        @test default_tortuosity.tortuosity_strength == 1.
        fixed_tortuosity = mr.fix(mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1.))],
            extracellular_fraction=0.2,
            tortuosity_strength=0.4,
        ))
        @test fixed_tortuosity.geometry.tortuosity_strength == 0.4
        for strength in (-0.1, NaN, Inf)
            @test_throws ArgumentError mr.fix(mr.LiminalGeometry(
                geometries=[(1., mr.Spheres(radius=1.))],
                extracellular_fraction=0.2,
                tortuosity_strength=strength,
            ))
        end
        # Strengths above one are valid if the morphology-dependent tensor is PSD.
        strong_tortuosity = mr.fix(mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1.))],
            extracellular_fraction=0.2,
            tortuosity_strength=2.,
        ))
        @test strong_tortuosity.geometry.tortuosity_strength == 2.
        @test_throws ArgumentError mr.fix(mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1.))],
            extracellular_fraction=0.2,
            tortuosity_strength=4.,
        ))
        for alpha in (0., 0.5, 1.)
            fixed = mr.fix(mr.LiminalGeometry(
                geometries=[(1., mr.Spheres(radius=1.))],
                extracellular_fraction=alpha,
            ))
            @test GI.tortuosity_tensor(fixed, GI.IsInside(())) ≈
                ((2 + alpha) / 3) * I atol=0.01
        end
        # A singular PSD tensor is supported, but a negative eigenvalue is not.
        LG = GI.PhysicalGeometries.LiminalGeometries
        _, transform, inverse_transform = LG._tortuosity_transforms([SA[1., 0., 0.]], 1.)
        @test transform == Diagonal([0., 1., 1.])
        @test inverse_transform == Diagonal([0., 1., 1.])
        @test_throws ArgumentError LG._tortuosity_transforms([SA[1., 0., 0.]], 1.1)
        fixed_liminal = mr.fix(mr.LiminalGeometry(
            geometries=[(1., mr.Spheres(radius=1.))],
            extracellular_fraction=0.2,
            tortuosity_strength=0.4,
        ))
        expected_tensor = I - 0.4 * (1 - 0.2) * fixed_liminal.geometry.normal_second_moment_tensor
        @test GI.tortuosity_tensor(
            fixed_liminal,
            GI.IsInside(()),
        ) ≈ expected_tensor
        @test GI.tortuosity_tensor(
            fixed_liminal,
            GI.IsInside((1,)),
        ) == I
        @test GI.tortuosity_tensor(
            mr.fix(mr.Spheres(radius=1.)),
            GI.IsInside(()),
        ) == I

        heterogeneous = GI.PhysicalGeometries.LiminalGeometries.FixedLiminalGeometry(
            [BaseObstructions.Sphere(1.), BaseObstructions.FiniteCylinder(
                SVector(0., 0., 0.), SVector(0., 0., 1.), 1.,
            )],
            [0.5, 0.5],
            0.,
        )
        @test eltype(heterogeneous.geometries) == Union{
            BaseObstructions.Sphere, BaseObstructions.FiniteCylinder,
        }
        @test GI.PhysicalGeometries.child_type(typeof(heterogeneous)) ===
            Core.apply_type(Transformations.Shift, 3, eltype(heterogeneous.geometries))
        @test GI.PhysicalGeometries.has_inside(typeof(heterogeneous))
        @test !GI.PhysicalGeometries.has_single_inside(typeof(heterogeneous))
        @test GI.PhysicalGeometries.inside_indices_eltype(typeof(heterogeneous)) ==
            Tuple{Tuple{Int, SVector{3, Float64}}}
        @test GI.PhysicalGeometries.intersection_type(typeof(heterogeneous)) ==
            Union{
                Tuple{Tuple{Int, SVector{3, Float64}}},
                Tuple{Tuple{Int, SVector{3, Float64}}, Int},
            }
    end

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
    @test_throws mr.BoundingBoxNotSupported BoundingBoxes.InternalBoundingBox(cylinder_rotation)
    collapsed_cylinder_box = BoundingBoxes.InternalBoundingBox(
        cylinder_rotation; no_deproject=true,
    )
    @test all(isfinite, BoundingBoxes.lower(collapsed_cylinder_box))
    @test all(BoundingBoxes.upper(collapsed_cylinder_box) .>
        BoundingBoxes.lower(collapsed_cylinder_box))
    no_deproject_positions, _ = GI.PhysicalGeometries.random_surface_positions(
        cylinder_rotation,
        Properties.GeometryLeafProperties(1.),
        collapsed_cylinder_box,
        1.;
        no_deproject=true,
    )
    @test all(iszero(position[3]) for position in no_deproject_positions)
    tilted_cylinder = Rotate(cylinder, [1 / sqrt(2) 0.; 0. 1.; 1 / sqrt(2) 0.])
    Random.seed!(1234)
    tilted_surface = GI.PhysicalGeometries.estimate_surface(
        tilted_cylinder; density=100., minimum_samples=1_000,
    )
    @test tilted_surface.area ≈ 2π rtol=0.1
    tilted_volume = GI.PhysicalGeometries.estimate_volume(
        cylinder_rotation; nsamples=100_000, rng=Random.MersenneTwister(1234),
    )
    @test tilted_volume ≈ π atol=0.05
    liminal_cylinder = GI.PhysicalGeometries.LiminalGeometries.FixedLiminalGeometry(
        [cylinder_rotation], [1.], 0.2,
    )
    @test !isempty(liminal_cylinder.outer_surface_sampling.positions)
    GI.PhysicalGeometries.LiminalGeometries.sample!(
        liminal_cylinder.outer_surface_sampling, liminal_cylinder, 100,
    )
    @test !isempty(liminal_cylinder.outer_surface_sampling.positions)
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
    repeat_distance = 2.0
    shifted_box_mesh = Mesh(
        [
            SVector(0.0, 0.0, 0.0),
            SVector(repeat_distance, 0.0, 0.0),
            SVector(0.0, repeat_distance, 0.0),
            SVector(repeat_distance, repeat_distance, 0.0),
            SVector(0.0, 0.0, repeat_distance),
            SVector(repeat_distance, 0.0, repeat_distance),
            SVector(0.0, repeat_distance, repeat_distance),
            SVector(repeat_distance, repeat_distance, repeat_distance),
        ],
        [
            SVector(1, 2, 3), SVector(2, 3, 4),
            SVector(5, 6, 7), SVector(6, 7, 8),
            SVector(1, 7, 5), SVector(1, 3, 7),
            SVector(2, 8, 6), SVector(2, 4, 8),
            SVector(1, 2, 6), SVector(1, 5, 6),
            SVector(3, 4, 8), SVector(3, 7, 8),
        ],
    )
    repeated_shifted_box = Repeat(
        shifted_box_mesh,
        fill(repeat_distance, 3),
    )
    @test inside_indices(repeated_shifted_box, SVector(1.5, 1.5, 1.5)) ==
        [(SVector(0, 0, 0),)]
    @test inside_indices(repeated_shifted_box, SVector(2.5, 1.5, 1.5)) ==
        [(SVector(1, 0, 0),)]
    @test inside_indices(repeated_shifted_box, SVector(-0.5, 1.5, 1.5)) ==
        [(SVector(-1, 0, 0),)]
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
    @test Set(mesh.indices[mesh.first_index_of_gap]) == Set((4, 2, 3))
    @test BoundingBoxes.lower(mesh.bounding_box) == SVector(0.0, 0.0, 0.0)
    @test BoundingBoxes.upper(mesh.bounding_box) == SVector(1.0, 1.0, 1.0)
    surface_density = Properties.GeometryLeafProperties(1.)
    _, surface_indices_without_gap = GI.PhysicalGeometries.random_surface_positions(
        mesh,
        surface_density,
        mesh.bounding_box,
        10_000;
        include_gap=false,
    )
    _, surface_indices_with_gap = GI.PhysicalGeometries.random_surface_positions(
        mesh,
        surface_density,
        mesh.bounding_box,
        10_000;
        include_gap=true,
    )
    @test all(index[1] < mesh.first_index_of_gap for index in surface_indices_without_gap)
    @test any(index[1] >= mesh.first_index_of_gap for index in surface_indices_with_gap)
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
    @test BoundingBoxes.lower(repeated_sphere.normalized_bounding_box) == SVector(-0.25, -0.25, -0.25)
    @test BoundingBoxes.upper(repeated_sphere.normalized_bounding_box) == SVector(0.25, 0.25, 0.25)
    @test inside_indices(repeated_sphere, SVector(4.5, 0.0, 0.0)) == [(SVector(1, 0, 0),)]
    @test inside_indices(repeated_sphere, SVector(2.0, 0.0, 0.0)) == Tuple{SVector{3, Int}}[]
    shifted_repeated_sphere = Repeat(
        Shift(BaseObstructions.Sphere(0.6), SVector(-1.0, 0.0, 0.0)),
        [2.0, 4.0, 4.0],
    )
    shifted_position = SVector(-0.5, 0.0, 0.0)
    @test inside_indices(shifted_repeated_sphere, shifted_position) == [(SVector(0, 0, 0),)]
    neighboring_position = SVector(0.5, 0.0, 0.0)
    @test Repeats.first_repeat(shifted_repeated_sphere, neighboring_position) == SVector(1, 0, 0)
    @test Repeats.last_repeat(shifted_repeated_sphere, neighboring_position) == SVector(1, 0, 0)
    @test inside_indices(shifted_repeated_sphere, neighboring_position) == [(SVector(1, 0, 0),)]
    @test_throws ArgumentError Repeat(BaseObstructions.Sphere(1.0), [0.0, 4.0, 4.0])
    @test_throws ArgumentError Repeat(BaseObstructions.Sphere(2.0), [1.0, 4.0, 4.0])
    @test_throws ArgumentError Repeat(
        Shift(BaseObstructions.Sphere(1.0), [3.0, 0.0, 0.0]),
        [2.0, 4.0, 4.0],
    )
        @test_throws mr.BoundingBoxNotSupported BoundingBoxes.InternalBoundingBox(repeated_sphere)
    repeated_hit = GI.PhysicalGeometries.find_intersection(
        repeated_sphere,
        SVector(-5.0, 0.0, 0.0),
        SVector(5.0, 0.0, 0.0),
    )
    @test repeated_hit[end] ≈ 0.4
    overlapping_repeat = Repeat(BaseObstructions.Sphere(1.5), [2.0, 4.0, 4.0])
    @test BoundingBoxes.lower(overlapping_repeat.normalized_bounding_box) == SVector(-0.75, -0.375, -0.375)
    @test BoundingBoxes.upper(overlapping_repeat.normalized_bounding_box) == SVector(0.75, 0.375, 0.375)
    @test inside_indices(overlapping_repeat, SVector(1.4, 0.0, 0.0)) == [
        (SVector(0, 0, 0),), (SVector(1, 0, 0),),
    ]
    overlapping_hit = GI.PhysicalGeometries.find_intersection(
        overlapping_repeat,
        SVector(1.4, 0.0, 0.0),
        SVector(1.9, 0.0, 0.0),
    )
    @test overlapping_hit[end] ≈ 0.2

    Random.seed!(1234)
    sampled_repeat = Repeat(
        GeometryVector([BaseObstructions.Sphere(0.1)]),
        [2.0, 4.0, 4.0],
    )
    sampled_box = BoundingBoxes.InternalBoundingBox(
        [1.1, 1.1, 1.1],
        [1.0, 0.0, 0.0],
    )
    sampled_positions, sampled_indices = GI.PhysicalGeometries.random_surface_positions(
        sampled_repeat,
        Properties.GeometryLeafProperties(1.0),
        sampled_box,
        100.0,
    )
    @test !isempty(sampled_positions)
    @test all(
        index[1][1] in (0, 1) && index[1][2:3] == SVector(0, 0)
        for index in sampled_indices
    )

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
    params = GI.PhysicalGeometries.IntersectionParams{3}(true, normal, false)
    @test shift isa Transformations.Transformation{3, 3, typeof(geometry_3d)}
    @test Transformations.to_child_coordinates(shift, position) == SVector(3.0, 3.0, 3.0)
    @test Transformations.from_child_coordinates(shift, SVector(3.0, 3.0, 3.0)) == position
    @test Transformations.to_child_coordinates_normal(shift, normal) == normal
    @test Transformations.from_child_coordinates(shift, params) === params

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
    @test_throws mr.BoundingBoxNotSupported BoundingBoxes.InternalBoundingBox(GeometryVector{3}(TestGeometry{3}[]))
    @test_throws mr.BoundingBoxNotSupported BoundingBoxes.InternalBoundingBox(GeometryTuple{3}(()))
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
    @test GI.PhysicalGeometries.from_child_coordinates(scale, params) === params
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
    projection_params = GI.PhysicalGeometries.IntersectionParams{2}(true, position_2d, false)
    @test projection isa Transformations.Transformation{3, 2, typeof(geometry_2d)}
    @test Transformations.to_child_coordinates(projection, SVector(1.0, 2.0, 3.0)) == SVector(1.0, 2.0)
    @test Transformations.from_child_coordinates(projection, SVector(1.0, 2.0)) == SVector(1.0, 2.0, 0.0)
    @test_throws ArgumentError GI.PhysicalGeometries.to_child_coordinates_normal(projection, SVector(1.0, 2.0, 3.0))
    projected_params = GI.PhysicalGeometries.from_child_coordinates(projection, projection_params)
    @test projected_params.inside
    @test projected_params.normal ≈ SVector(1.0, 2.0, 0.0) / sqrt(5)
    @test !projected_params.hit_gap
    @test_throws mr.BoundingBoxNotSupported BoundingBoxes.InternalBoundingBox(projection)

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
    box_1d = BoundingBoxes.InternalBoundingBox{1}(1.0)
    @test BoundingBoxes.could_intersect(box_1d, SVector(-2.0), SVector(2.0))
    @test !BoundingBoxes.could_intersect(box_1d, SVector(2.0), SVector(3.0))
    box_2d = BoundingBoxes.InternalBoundingBox{2}(1.0)
    @test BoundingBoxes.could_intersect(box_2d, SVector(-2.0, 0.0), SVector(2.0, 0.0))
    @test !BoundingBoxes.could_intersect(box_2d, SVector(2.0, 2.0), SVector(3.0, 3.0))
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

@testset "Surface distance" begin
    position = SVector(0.0, 0.0, 0.0)
    @test isinf(GI.distance_to_surface(mr.fix(()), position))

    wall = mr.fix(mr.Walls(position=0.0))
    @test GI.distance_to_surface(wall, SVector(1.0, 0.0, 0.0)) == 1.0
    @test GI.distance_to_surface(wall, SVector(-2.0, 0.0, 0.0)) == 2.0

    sphere = mr.fix(mr.Spheres(radius=1.0))
    @test GI.distance_to_surface(sphere, position) == 1.0
    @test GI.distance_to_surface(sphere, SVector(1.5, 0.0, 0.0)) == 0.5

    triangle = BaseObstructions.FullTriangle(
        SVector(0.0, 0.0, 0.0),
        SVector(1.0, 0.0, 0.0),
        SVector(0.0, 1.0, 0.0),
    )
    @test GI.distance_to_surface(triangle, SVector(0.25, 0.25, 1.0)) == 1.0

    repeated_wall = mr.fix(mr.Walls(position=0.0, repeats=2.0))
    @test GI.distance_to_surface(repeated_wall, SVector(2.1, 0.0, 0.0)) ≈ 0.1
    @test GI.distance_to_surface(repeated_wall, SVector(102.1, 0.0, 0.0)) ≈ 0.1

    combined = mr.fix([
        mr.Spheres(radius=1.0),
        mr.Spheres(position=[5.0, 0.0, 0.0], radius=1.0),
    ])
    @test GI.distance_to_surface(combined, SVector(1.5, 0.0, 0.0)) == 0.5
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
    @test GI.susceptibility_off_resonance(walls, susceptibility_position) == 0.0
    @test GI.off_resonance_gradient(walls, 1.0) == 0.0

    z_cylinder = mr.Cylinders(
        radius=1.0,
        g_ratio=0.8,
        susceptibility_iso=1.0,
        susceptibility_aniso=0.0,
        grid_resolution=Inf,
    )
    fixed_z_cylinder = mr.fix(z_cylinder)
    @test fixed_z_cylinder.susceptibility == ()
    @test GI.susceptibility_off_resonance(fixed_z_cylinder, susceptibility_position) == 0.0

    cylinder = mr.Cylinders(
        radius=1.0,
        rotation=:x,
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
    @test GI.susceptibility_off_resonance(fixed_cylinder, susceptibility_position) ==
        GI.Susceptibility.susceptibility_off_resonance(
            fixed_cylinder.susceptibility,
            susceptibility_position,
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

    sphere = mr.Spheres(
        radius=1.0,
        susceptibility=1.0,
        grid_resolution=Inf,
    )
    liminal = mr.LiminalGeometry(
        geometries=[(1.0, sphere), (1.0, cylinder)],
        extracellular_fraction=0.2,
    )
    fixed_liminal = mr.fix(liminal)
    @test fixed_liminal.susceptibility isa GI.Susceptibility.Liminal.LiminalSusceptibility
    @test GI.susceptibility_off_resonance(
        fixed_liminal, susceptibility_position,
    ) == 0.0
    offset = SVector(1.0, 2.0, 3.0)
    sphere_hit = GI.Intersection(
        0.0, ((1, offset),), (), zero(SVector{3, Float64}), true, false,
    )
    cylinder_hit = GI.Intersection(
        0.0, ((2, offset),), (), zero(SVector{3, Float64}), true, false,
    )
    sphere_inside = GI.IsInside([((1, offset),)])
    cylinder_inside = GI.IsInside([((2, offset),)])
    sphere_field = GI.susceptibility_off_resonance(
        mr.fix(sphere), susceptibility_position - offset,
    )
    sphere_surface_field = GI.susceptibility_off_resonance(
        mr.fix(sphere), susceptibility_position - offset, sphere_hit,
    )
    cylinder_field = GI.susceptibility_off_resonance(
        fixed_cylinder, susceptibility_position - offset,
    )
    @test GI.susceptibility_off_resonance(
        fixed_liminal, susceptibility_position, sphere_hit, sphere_inside,
    ) == sphere_field
    @test GI.susceptibility_off_resonance(
        fixed_liminal, susceptibility_position, nothing, sphere_inside,
    ) == sphere_field
    @test GI.susceptibility_off_resonance(
        fixed_liminal, susceptibility_position, nothing, GI.IsInside(()),
    ) == 0.0
    @test GI.susceptibility_off_resonance(
        fixed_liminal, susceptibility_position, sphere_hit, sphere_inside,
    ) == sphere_surface_field
    @test GI.susceptibility_off_resonance(
        fixed_liminal, susceptibility_position, cylinder_hit, cylinder_inside,
    ) == cylinder_field
    @test GI.off_resonance_gradient(fixed_liminal, 1.0) == max(
        GI.off_resonance_gradient(mr.fix(sphere), 1.0),
        GI.off_resonance_gradient(fixed_cylinder, 1.0),
    )
end
