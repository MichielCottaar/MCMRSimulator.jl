"""
Types:
- [`FixedGeometry`](@ref)
- [`FixedObstructionGroup`](@ref)

Methods:
- [`has_inside`](@ref)
- [`isinside`](@ref)
- [`detect_intersection`](@ref)
- [`size_scale`](@ref)
"""
module FixedObstructionGroups

import Random: rand
import LinearAlgebra: inv, transpose, norm, cross, ⋅, I
import StaticArrays: SVector, SMatrix
import ..Obstructions:
    FixedObstruction, has_inside, isinside, obstruction_type, random_surface_positions,
    IndexTriangle, FullTriangle, size_scale, Shift, Wall, curvature
import ..BoundingBoxes: BoundingBox, could_intersect
import ..Intersections: Intersection, empty_intersection
import ..Reflections: Reflection, empty_reflection
import ..HitGrids: HitGrid, detect_intersection_grid


"""
Collection of multiple base [`FixedObstruction`](@ref) objects.

This is the main internal representation of a group of identical [`FixedObstruction`](@ref) objects.

Properties:
- `obstructions`: vector of the actual [`FixedObstruction`](@ref) objects.
- `parent_index`: Index of this group within the larger [`FixedGeometry`](@ref).
- `original_index`: Index of this group within the original user-provided geometry.
- `rotation`: rotation from global 3-dimensional space to the 1, 2, or 3-dimensional space of the obstructions.
- `inv_rotation`: inverse of the rotation above
- `grid`: [`Grid`](@ref) object on which the obstruction intersections have been precomputed. This speeds up the detection of intersections.
- `bounding_boxes`: vector of [`BoundingBox`](@ref) objects for each obstruction. These are used to predect whether a spin could intersect with the obstruction.
- `volume`: R1, R2, and off-resonance properties of the spins inside the obstructions.
- `surface`: R1, R2, off-resonance, surface_density and dwell_time properties of particles stuck to the surface. Also, contains the permeability and surface relaxivity to process collsions.
- `vertices`: vector of vertices (only used for a mesh).
"""
struct FixedObstructionGroup{
    N, 
    R <: Union{Nothing, SVector{N, Float64}}, 
    O <: FixedObstruction{N},
    G <: HitGrid{N, O},
    V <: NamedTuple{(:R1, :R2, :off_resonance)},
    S <: NamedTuple{(:R1, :R2, :off_resonance, :permeability, :surface_density, :dwell_time, :surface_relaxivity)},
    A <: Tuple,
    K
    }
    repeats :: R
    parent_index :: Int
    original_index :: Int

    # rotations
    rotation :: SMatrix{3, N, Float64, K}
    inv_rotation :: SMatrix{N, 3, Float64, K}

    # obstruction indices
    hit_grid :: G

    # MRI and collision properties
    # single value, or one value per obstruction stored in NamedTuple
    volume :: V
    surface :: S

    size_scale :: Float64

    # Additional arguments that should be passed around for the `obstructions` to work
    args :: A
    function FixedObstructionGroup(obstructions, repeats, parent_index, original_index, rotation, grid, volume, surface, size_scale, args...)
        new{
            size(rotation, 2), typeof(repeats), eltype(obstructions), typeof(grid),
            typeof(volume), typeof(surface), typeof(args), 3 * size(rotation, 2)
        }(
            obstructions, repeats, 
            parent_index, original_index, 
            rotation, transpose(rotation), 
            grid, volume, surface,
            size_scale, args,
        )
    end
end

"""
    FixedGeometry([obstruction_groups...])

A collection of [`FixedObstructionGroup`](@ref) objects each reperesenting part of the geometry.
"""
const FixedGeometry{N} = NTuple{N, FixedObstructionGroup}

repeating(::FixedObstructionGroup{N, R}) where {N, R} = R <: SVector{N, Float64}
repeating(::Type{<:FixedObstructionGroup{N, R}}) where {N, R} = R <: SVector{N, Float64}

obstruction_type(::Type{<:FixedObstructionGroup{N, R, O}}) where {N, R, O} = obstruction_type(O)

rotate_from_global(g::FixedObstructionGroup, pos::SVector{3}) = g.inv_rotation * pos
rotate_to_global(g::FixedObstructionGroup{N}, pos::SVector{N}) where {N} = g.rotation * pos


curvature(g::FixedMesh) = length(g.obstructions) == 1 ? 0. : curvature(g.obstructions, g.vertices)


size_scale(g::FixedObstructionGroup) = g.size_scale
size_scale(g::FixedGeometry) = minimum(size_scale.(g))
size_scale(g::FixedGeometry{0}) = Inf

function Base.show(io::IO, geom::FixedObstructionGroup)
    print(io, length(geom.obstructions), " ")
    if repeating(geom)
        print(io, "repeating ")
    end
    print(io, String(nameof(obstruction_type(typeof(geom)))) * " objects")
end

"""
    BoundingBox(obstruction_group)

Finds the bounding box containing all the obstructions in the [`FixedObstructionGroup`](@ref) ignoring any repeats or rotation.
"""
BoundingBox(group::FixedObstructionGroup) = BoundingBox(group.obstructions, group.vertices)
BoundingBox(obstructions::Vector{<:FixedObstruction}, vertices::Vector{SVector{3, Float64}}) = BoundingBox(BoundingBox.(obstructions))
BoundingBox(obstructions::Vector{IndexTriangle}, vertices::Vector{SVector{3, Float64}}) = BoundingBox([BoundingBox(o, vertices) for o in obstructions])

"""
    isinside(obstruction_group, position[, stuck_to])

Returns a vector of indices with all the obstructions in [`FixedObstructionGroup`](@ref) containing the `position` (in order).
For obstructions with only a single inside, will return an empty vector ("[]") if the particle is outside and a "[0]" if inside.
"""
function isinside(g::FixedObstructionGroup, pos::SVector{3}, stuck_to::Reflection=empty_reflection)
    isinside(g, pos, stuck_to.geometry_index == g.parent_index ? stuck_to.obstruction_index : 0, stuck_to.inside)
end

function isinside(g::FixedObstructionGroup, pos::SVector{3}, stuck_to, inside::Bool)
    if ~has_inside(typeof(g))
        return Int[]
    end
    rotated = rotate_from_global(g, pos)

    if repeating(g)
        repeats = g.grid.size
        half_repeats = repeats/2
        voxel = @. Int(div(rotated, repeats, RoundNearest))
        normed = rotated .- voxel .* repeats
    else
        normed = rotated
    end
    indices = Int[]
    for (index, shift) in get_indices(g.grid, normed)
        if stuck_to == index
            if inside
                push!(indices, index)
            end
            continue
        end
        if iszero(shift)
            shifted = normed
        else
            shifted = normed .- g.grid.shifts[shift]
        end
        if isinside(g.obstructions[index], shifted)
            push!(indices, index)
        end
    end
    return indices
end


"""
    isinside(geometry, position[, stuck_reflection])

Returns a vector of pairs of indices with all the obstructions in [`FixedGeometry`](@ref) containing the `position` (in order).
The first index indicates the index of the [`FixedObstructionGroup`](@ref) within the `geometry`.
The second index is the index of the specific [`FixedObstruction`](@ref) within the obstruction group.
"""
function isinside(gs::FixedGeometry, pos::SVector{3}, stuck_to::Reflection=empty_reflection)
    indices = Tuple{Int, Int}[]
    for g in gs
        for inner_index in isinside(g, pos, stuck_to)
            push!(indices, (g.parent_index, inner_index))
        end
    end
    return indices
end


"""
    detect_intersection(obstruction_group/geometry, start, dest[, previous_intersection])

Find the closest intersection between the line from `start` to `dest` and an obstruction in the `geometry(ies)`.
"""
function detect_intersection(g::FixedObstructionGroup{N}, start::SVector{3}, dest::SVector{3}, previous_hit::Tuple{Int, Int, Bool}=(0, 0, false)) where {N}
    rotated_start = rotate_from_global(g, start)
    rotated_dest = rotate_from_global(g, dest)

    if previous_hit[1] != g.parent_index
        previous_index = 0
        prev_inside = false
    else
        previous_index = previous_hit[2]
        prev_inside = previous_hit[3]
    end
    if repeating(g)
        (index, intersection) = detect_intersection_repeating(g, rotated_start, rotated_dest, previous_index, prev_inside)
    else
        (index, intersection) = detect_intersection_grid(g.grid, rotated_start, rotated_dest, previous_index, prev_inside)
    end
    if intersection.distance > 1.
        return empty_intersection
    end
    return Intersection(
        intersection.distance,
        rotate_to_global(g, intersection.normal),
        intersection.inside,
        g.parent_index,
        index
    )
end

function detect_intersection_repeating(g::FixedObstructionGroup{N}, start::SVector{N}, dest::SVector{N}, prev_index::Int, prev_inside::Bool) where {N}
    grid_start = start ./ g.repeats
    grid_dest = dest ./ g.repeats
    voxel_f = round.(grid_start)
    if all(voxel_f .== round.(grid_dest))
        voxel_f2 = voxel_f .* g.repeats
        return detect_intersection_grid(g.grid, start .- voxel_f2, dest .- voxel_f2, prev_index, prev_inside)
    end
    found_index = zero(Int32)
    found_intersection = empty_obstruction_intersections[N]
    for (voxel, _, _, t2, _) in ray_grid_intersections(grid_start .+ 0.5, grid_dest .+ 0.5)
        scaled_voxel = voxel .* g.repeats
        (index, intersection) = detect_intersection_grid(g.grid, start .- scaled_voxel, dest .- scaled_voxel, prev_index, prev_inside)
        if found_intersection.distance > intersection.distance
            found_intersection = intersection
            found_index = index
        end
        if found_intersection.distance < t2
            return (found_index, found_intersection)
        end
    end
    return (found_index, found_intersection)
end


function detect_intersection(geometries::FixedGeometry, start::SVector{3}, dest::SVector{3}, previous_intersection=(0, 0, false))
    intersection = empty_intersection
    for geometry in geometries
        new_intersection = detect_intersection(geometry, start, dest, previous_intersection)
        if new_intersection.distance < intersection.distance
            intersection = new_intersection
        end
    end
    return intersection
end


"""
    random_surface_positions(group/geometry, bounding_box, volume_density)

Randomly draws positions on the surface within the [`BoundingBox`](@ref).
The density of points will be equal to `surface_density` * `volume_density`,
where `surface_density` is either set by the `group` itself or by `default_surface_density` (which should come from [`GlobalProperties`](@ref)).

For each drawn position will return a tuple with:
- position: 3D in um
- surface normal: 3D, normalised (always pointing inwards)
- geometry_index: index of the group
- obstruction_index: index of the obstruction within the group that the position is on the surface of
"""
function random_surface_positions(group::FixedObstructionGroup{N}, bb::BoundingBox{3}, volume_density::Number) where {N}
    local_surface_density = group.surface.surface_density
    if local_surface_density isa Vector
        surface_density = [isnothing(sd) ? default_surface_density : sd for sd in local_surface_density]
    else
        surface_density = fill(local_surface_density, length(group.obstructions))
    end

    if N == 1
        normal = group.rotation[:, 1]
        e1 = cross(SVector{3}([1., 1, 0]), normal)
        if all(iszero.(e1))
            e1 = cross(SVector{3}([1., 0, 0]), normal)
        end
        e1 = e1 ./ norm(e1)
        e2 = cross(normal, e1)
        bb_area = (abs.(e1) ⋅ (bb.upper - bb.lower)) * (abs.(e2) ⋅ (bb.upper - bb.lower))
    elseif N == 2
        n1 = group.rotation[:, 1]
        n2 = group.rotation[:, 2]
        edge = cross(n1, n2)
        bb_area = abs.(edge) ⋅ (bb.upper - bb.lower)
    else
        bb_area = 1
    end
    normed_surface_density = surface_density .* (bb_area * volume_density)

    get_random_pos(o, d) = random_surface_positions(o, group.vertices, d)

    if repeating(group)
        repeats = group.repeats
        normals = [group.rotation[:, i] for i in 1:N]
        nrepeats_lower = div.([bb.lower ⋅ abs.(n) for n in normals], group.grid.size, RoundDown) .- 5
        nrepeats_upper = div.([bb.upper ⋅ abs.(n) for n in normals], group.grid.size, RoundUp) .+ 5

        total_repeats = prod(nrepeats_upper .- nrepeats_lower .+ 1)
        sub_draws = get_random_pos.(group.obstructions, normed_surface_density .* total_repeats)
        base_positions = [p for s in sub_draws for p in s[1]]
        normals = [p for s in sub_draws for p in s[2]]
        shifts = SVector{N, Float64}.(zip(rand.(UnitRange.(nrepeats_lower, nrepeats_upper), length(base_positions))...))
        positions = base_positions .+ [s .* repeats for s in shifts]
    else
        sub_draws = get_random_pos.(group.obstructions, normed_surface_density)
        positions = [p for s in sub_draws for p in s[1]]
        normals = [p for s in sub_draws for p in s[2]]
    end

    global_normals = [rotate_to_global(group, n) for n in normals]
    unrotated_positions = [rotate_to_global(group, p) for p in positions]
    if N == 3
        global_positions = unrotated_positions
    elseif N == 2
        l1 = abs.(edge) ⋅ bb.lower
        l2 = abs.(edge) ⋅ bb.upper
        global_positions = [(rand() * (l2 - l1) .+ l1) .* edge .+ p for p in unrotated_positions]
    else
        l1l = abs.(e1) ⋅ bb.lower
        l1u = abs.(e1) ⋅ bb.upper
        l2l = abs.(e2) ⋅ bb.lower
        l2u = abs.(e2) ⋅ bb.upper
        global_positions = [(rand() * (l1u - l1l) .+ l1l) .* e1 .+ (rand() * (l2u - l2l) .+ l2l) .* e2 .+ p for p in unrotated_positions]
    end
    geometry_index = fill(group.parent_index, length(positions))
    obstruction_index = vcat([fill(index, length(arr[1])) for (index, arr) in enumerate(sub_draws)]...)

    return filter(
        t->all(t[1] .>= bb.lower) && all(t[1] .<= bb.upper),
        collect(zip(global_positions, global_normals, geometry_index, obstruction_index))
    )
end

function random_surface_positions(geometry::FixedGeometry, bb::BoundingBox{3}, volume_density::Number)
    vcat([random_surface_positions(g, bb, volume_density) for g in geometry]...)
end

function random_surface_positions(geometry::FixedGeometry{0}, bb::BoundingBox{3}, volume_density::Number)
    return Tuple{SVector{3, Float64}, SVector{3, Float64}, Int{}, Int{}}[]
end

end