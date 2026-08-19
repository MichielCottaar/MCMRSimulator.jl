"""Physical geometries that consist of groups of other physical geometries."""
module Groups

import StaticArrays: SVector
import ...Indices: ObstructionIndex, add_index
import ...InternalBoundingBoxes
import ..GridDispatch: IntersectionGrid, GridIterator, detect_intersection_grid
import ..Intersections: remove_expected_index
import ..PhysicalGeometries: PhysicalGeometry, Intersection, find_intersection, has_inside, has_single_inside, isinside_single, inside_indices, InternalBoundingBox
import ..PhysicalGeometries: intersection_index_length, inside_index_length, all_equal_inside_depth
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh
import ...Properties: GeometryProperties, GeometryLeafProperties, GeometryVectorProperties, GeometryTupleProperties

abstract type GroupGeometry{N, P} <: PhysicalGeometry{N} end
abstract type GeometryVectorLike{N, P<:PhysicalGeometry{N}} <: GroupGeometry{N, P} end

group_geometries(group::GeometryVectorLike) = group.geometries


function find_intersection(group::GroupGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, previous_hit=nothing)
    current = nothing
    for (index, candidate, dist_all_checked) in intersection_candidates(group, start, dest)
        if !isnothing(current) && current[end] < dist_all_checked
            return current
        end
        candidate_previous_hit = isnothing(previous_hit) || previous_hit[1] != index ? nothing : previous_hit[2:end]
        intersect = find_intersection(candidate, start, dest, candidate_previous_hit)
        if isnothing(intersect)
            continue
        end
        if isnothing(current) || intersect[end] < current[end]
            current = (index, intersect...)
        end
    end
    return current
end

_is_fully_concrete(::Type{T}) where {T} =
    isconcretetype(T) && all(parameter -> !(parameter isa Type) || _is_fully_concrete(parameter), T.parameters)

struct GeometryVector{N, P<:PhysicalGeometry{N}} <: GeometryVectorLike{N, P}
    geometries::Vector{P}

    function GeometryVector{N, P}(geometries::Vector{P}) where {N, P<:PhysicalGeometry{N}}
        _is_fully_concrete(P) || throw(ArgumentError("GeometryVector element types must be fully concrete"))
        new{N, P}(geometries)
    end
end

struct GeometryVectorBoundingBox{N, P<:PhysicalGeometry{N}} <: GeometryVectorLike{N, P}
    geometries::Vector{P}
    bounding_boxes::Vector{InternalBoundingBoxes.InternalBoundingBox{N}}

    function GeometryVectorBoundingBox{N, P}(
        geometries::Vector{P},
        bounding_boxes::Vector{InternalBoundingBoxes.InternalBoundingBox{N}},
    ) where {N, P<:PhysicalGeometry{N}}
        length(geometries) == length(bounding_boxes) ||
            throw(ArgumentError("geometries and bounding_boxes must have the same length"))
        _is_fully_concrete(P) || throw(ArgumentError("GeometryVector element types must be fully concrete"))
        new{N, P}(geometries, bounding_boxes)
    end
end

struct GeometryVectorGrid{N, P<:PhysicalGeometry{N}} <: GeometryVectorLike{N, P}
    geometries::Vector{P}
    grid::IntersectionGrid{N}

    function GeometryVectorGrid{N, P}(
        geometries::Vector{P},
        grid::IntersectionGrid{N},
    ) where {N, P<:PhysicalGeometry{N}}
        _is_fully_concrete(P) || throw(ArgumentError("GeometryVector element types must be fully concrete"))
        new{N, P}(geometries, grid)
    end
end

struct GeometryTuple{N, P<:Tuple{Vararg{PhysicalGeometry{N}}}} <: GroupGeometry{N, P}
    geometries::P
end

group_geometries(group::GeometryTuple) = group.geometries

intersection_candidates(group::Union{GeometryVector, GeometryTuple}, start, destination) =
    ((index, child, 0.0) for (index, child) in enumerate(group_geometries(group)))

intersection_candidates(group::GeometryVectorBoundingBox, start, destination) =
    (
        (index, group_geometries(group)[index], 0.0)
        for index in eachindex(group_geometries(group))
        if _could_intersect(group, index, start, destination)
    )

intersection_candidates(group::GeometryVectorGrid, start, destination) =
    (
        (child_index, group_geometries(group)[child_index], dist_all_checked)
        for (child_index, dist_all_checked) in GridIterator(group.grid, start, destination)
    )

has_inside(::Type{<:GeometryVectorLike{N, P}}) where {N, P} = has_inside(P)
has_inside(::Type{<:GeometryTuple{N, P}}) where {N, P} = any(has_inside, P.parameters)
has_single_inside(::Type{<:GroupGeometry}) = false

_tuple_intersection_index_length(::Type{P}) where {P<:Tuple} =
    isempty(P.parameters) ? 0 : intersection_index_length(P.parameters[1]) + 1
_tuple_inside_index_length(::Type{P}) where {P<:Tuple} =
    isempty(P.parameters) ? 0 : inside_index_length(P.parameters[1]) + 1

intersection_index_length(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    _tuple_intersection_index_length(P)
inside_index_length(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    _tuple_inside_index_length(P)
intersection_index_length(::Type{<:GeometryVectorLike{N, P}}) where {N, P} =
    intersection_index_length(P) + 1
inside_index_length(::Type{<:GeometryVectorLike{N, P}}) where {N, P} =
    inside_index_length(P) + 1
all_equal_inside_depth(::Type{<:GeometryVectorLike{N, P}}) where {N, P} =
    all_equal_inside_depth(P)

function all_equal_inside_depth(::Type{<:GeometryTuple{N, P}}) where {N, P}
    child_types = P.parameters
    isempty(child_types) && return true
    all(all_equal_inside_depth, child_types) || return false
    first_depth = inside_index_length(child_types[1])
    all(inside_index_length(child_type) == first_depth for child_type in child_types)
end

@generated function _empty_inside_indices(::Type{G}) where {G}
    all_equal_inside_depth(G) || return :(ObstructionIndex[])
    depth = inside_index_length(G)
    :(ObstructionIndex{$depth}[])
end

function InternalBoundingBox(geometry::GroupGeometry{N, P}) where {N, P}
    isempty(geometry) && throw(ArgumentError("cannot construct a bounding box for an empty geometry group"))
    boxes = [InternalBoundingBoxes.InternalBoundingBox(child) for child in geometry]
    lower_bound = reduce((a, b) -> min.(a, b), (InternalBoundingBoxes.lower(box) for box in boxes))
    upper_bound = reduce((a, b) -> max.(a, b), (InternalBoundingBoxes.upper(box) for box in boxes))
    InternalBoundingBoxes.InternalBoundingBox(
        (upper_bound - lower_bound) / 2,
        (upper_bound + lower_bound) / 2,
    )
end

InternalBoundingBox(geometry::GeometryVectorGrid) = geometry.grid.bounding_box

function _inside_indices(
    geometry::GeometryVectorLike{N, P},
    child_indices,
    position::SVector{N, Float64},
    intersection::Intersection{3, M}=Intersection{3, inside_index_length(typeof(geometry))}(),
) where {N, P, M}
    has_inside(P) || return _empty_inside_indices(typeof(geometry))
    geometry_type = typeof(geometry)
    indices = _empty_inside_indices(geometry_type)
    if has_single_inside(P)
        for child_index in child_indices
            child_intersection = remove_expected_index(intersection, child_index)
            isinside_single(group_geometries(geometry)[child_index], position, child_intersection) || continue
            push!(indices, add_index(ObstructionIndex(), child_index))
        end
    else
        for child_index in child_indices
            child_intersection = remove_expected_index(intersection, child_index)
            for obstruction_index in inside_indices(group_geometries(geometry)[child_index], position, child_intersection)
                push!(indices, add_index(obstruction_index, child_index))
            end
        end
    end
    indices
end

function inside_indices_for_any_type(geometry, position, intersection)
    geometry_type = typeof(geometry)
    has_inside(geometry_type) || return _empty_inside_indices(geometry_type)
    if has_single_inside(geometry_type)
        index_type = inside_index_length(geometry_type)
        return isinside_single(geometry, position, intersection) ?
            [ObstructionIndex{index_type}()] : ObstructionIndex{index_type}[]
    end
    return inside_indices(geometry, position, intersection)
end

function inside_indices(
    geometry::GeometryVectorLike{N},
    position::SVector{N, Float64},
    intersection::Intersection{3, M}=Intersection{3, inside_index_length(typeof(geometry))}(),
) where {N, M}
    child_indices = if geometry isa GeometryVectorBoundingBox
        (child_index for child_index in eachindex(group_geometries(geometry))
            if InternalBoundingBoxes.isinside(geometry.bounding_boxes[child_index], position))
    else
        eachindex(group_geometries(geometry))
    end
    _inside_indices(geometry, child_indices, position, intersection)
end

function _grid_coordinate(
    geometry::GeometryVectorGrid{N},
    position::SVector{N, Float64},
) where {N}
    coordinate = Int.(floor.((position - InternalBoundingBoxes.lower(geometry.grid.grid_bounding_box)) .* geometry.grid.inv_resolution)) .+ 1
    any(coordinate .< 1) || any(coordinate .> size(geometry.grid.indices)) ? nothing : SVector{N, Int}(coordinate)
end

function inside_indices(
    geometry::GeometryVectorGrid{N},
    position::SVector{N, Float64},
    intersection::Intersection{3, M}=Intersection{3, inside_index_length(typeof(geometry))}(),
) where {N, M}
    coordinate = _grid_coordinate(geometry, position)
    isnothing(coordinate) && return _empty_inside_indices(typeof(geometry))
    _inside_indices(geometry, _grid_candidates(geometry, coordinate), position, intersection)
end

function _grid_candidates(geometry::GeometryVectorGrid, coordinate)
    (any(coordinate .< 1) || any(coordinate .> size(geometry.grid.indices))) && return Int[]
    geometry.grid.indices[coordinate...]
end

function inside_indices(
    geometry::GeometryTuple{N},
    position::SVector{N, Float64},
    intersection::Intersection{3, M}=Intersection{3, inside_index_length(typeof(geometry))}(),
) where {N, M}
    geometry_type = typeof(geometry)
    indices = _empty_inside_indices(geometry_type)
    for (child_index, child) in enumerate(geometry)
        child_intersection = remove_expected_index(intersection, child_index)
        append!(indices, add_index.(inside_indices_for_any_type(child, position, child_intersection), child_index))
    end
    indices
end

GeometryVector{N}(geometries::Vector{P}) where {N, P<:PhysicalGeometry{N}} =
    GeometryVector{N, P}(geometries)

function GeometryVector(
    geometries::Vector{P};
    bounding_box::Bool=false,
    grid::Bool=false,
    grid_resolution=nothing,
) where {N, P<:PhysicalGeometry{N}}
    bounding_box && grid && throw(ArgumentError("bounding_box and grid cannot both be enabled"))
    grid && return GeometryVectorGrid(geometries; grid_resolution)
    bounding_box ? GeometryVectorBoundingBox(geometries) : GeometryVector{N, P}(geometries)
end

function GeometryVectorGrid(
    geometries::Vector{P};
    grid_resolution=nothing,
) where {N, P<:PhysicalGeometry{N}}
    isempty(geometries) && throw(ArgumentError("cannot construct a grid for an empty geometry vector"))
    boxes = InternalBoundingBox.(geometries)
    grid = IntersectionGrid(boxes; resolution=grid_resolution)
    GeometryVectorGrid{N, P}(geometries, grid)
end

GeometryVectorBoundingBox(geometries::Vector{P}) where {N, P<:PhysicalGeometry{N}} =
    GeometryVectorBoundingBox{N, P}(geometries, InternalBoundingBox.(geometries))

GeometryVectorBoundingBox(
    geometries::Vector{P},
    bounding_boxes::Vector{InternalBoundingBoxes.InternalBoundingBox{N}},
) where {N, P<:PhysicalGeometry{N}} =
    GeometryVectorBoundingBox{N, P}(geometries, bounding_boxes)

GeometryTuple{N}(geometries::Tuple{Vararg{PhysicalGeometry{N}}}) where {N} =
    GeometryTuple{N, typeof(geometries)}(geometries)

function _could_intersect(
    ::GroupGeometry,
    ::Int,
    ::SVector,
    ::SVector,
)
    true
end

function _could_intersect(
    geometry::GeometryVectorBoundingBox{N},
    child_index::Int,
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
 ) where {N}
    box = geometry.bounding_boxes[child_index]
    InternalBoundingBoxes.could_intersect(box, start, destination) &&
        InternalBoundingBoxes.does_intersect(box, start, destination)
end

function detect_intersection(
    geometry::GroupGeometry{N, P},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3, M}=Intersection{3, intersection_index_length(typeof(geometry))}(),
) where {N, P, M}
    if !Base.isempty(previous_hit)
        indices = previous_hit.obstruction_index.indices
        isempty(indices) && throw(ArgumentError("a non-empty previous hit must have an obstruction index"))
        1 <= indices[1] <= length(geometry) ||
            throw(ArgumentError("previous-hit obstruction index is not a child of this group"))
    end

    closest = Intersection{N, intersection_index_length(typeof(geometry))}()
    for (child_index, child) in enumerate(geometry)
        _could_intersect(geometry, child_index, start, destination) || continue
        child_previous_hit = remove_expected_index(previous_hit, child_index)
        intersection = detect_intersection(child, start, destination, child_previous_hit)
        if intersection.distance < closest.distance
            closest = add_index(intersection, child_index)
        end
    end
    return closest
end

function detect_intersection(
    geometry::GeometryVectorGrid{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3, M}=Intersection{3, intersection_index_length(typeof(geometry))}(),
) where {N, M}
    return detect_intersection_grid(
        geometry.grid,
        start,
        destination,
        previous_hit,
        Intersection{N, intersection_index_length(typeof(geometry))}(),
    ) do child_index, previous
        child_previous_hit = remove_expected_index(previous, child_index)
        intersection = detect_intersection(group_geometries(geometry)[child_index], start, destination, child_previous_hit)
        add_index(intersection, child_index)
    end
end

Base.size(geometry::GeometryVectorLike) = size(group_geometries(geometry))
Base.axes(geometry::GeometryVectorLike) = axes(group_geometries(geometry))
Base.length(geometry::GeometryVectorLike) = length(group_geometries(geometry))
Base.firstindex(geometry::GeometryVectorLike) = firstindex(group_geometries(geometry))
Base.lastindex(geometry::GeometryVectorLike) = lastindex(group_geometries(geometry))
Base.getindex(geometry::GeometryVectorLike, index...) = getindex(group_geometries(geometry), index...)
Base.setindex!(geometry::GeometryVectorLike, value, index...) = setindex!(group_geometries(geometry), value, index...)
Base.iterate(geometry::GeometryVectorLike, state...) = iterate(group_geometries(geometry), state...)
Base.eltype(::Type{<:GeometryVectorLike{N, P}}) where {N, P} = P

Base.length(geometry::GeometryTuple) = length(group_geometries(geometry))
Base.firstindex(geometry::GeometryTuple) = firstindex(group_geometries(geometry))
Base.lastindex(geometry::GeometryTuple) = lastindex(group_geometries(geometry))
Base.getindex(geometry::GeometryTuple, index...) = getindex(group_geometries(geometry), index...)
Base.iterate(geometry::GeometryTuple, state...) = iterate(group_geometries(geometry), state...)
Base.eltype(::Type{GeometryTuple{N, P}}) where {N, P} = eltype(P)
Base.Tuple(geometry::GeometryTuple) = group_geometries(geometry)

_density_child(density::GeometryLeafProperties, ::Int) = density
_density_child(density::GeometryVectorProperties, index::Int) = density.properties[index]
_density_child(density::GeometryTupleProperties, index::Int) = density.properties[index]

function _combine(draws, ::Val{N}) where {N}
    draws = collect(draws)
    positions = reduce(vcat, (draw[1] for draw in draws); init=SVector{N, Float64}[])
    intersections = reduce(vcat, (draw[2] for draw in draws); init=Intersection{N}[])
    positions, intersections
end

function random_surface_positions(geometry::Union{GeometryVectorLike{N}, GeometryTuple{N}}, density::GeometryProperties,
    bounding_box::InternalBoundingBox{N}, scale_density) where {N}
    _combine((let
        values = random_surface_positions(child, _density_child(density, index), bounding_box, scale_density)
        (values[1], [add_index(hit, index) for hit in values[2]])
    end for (index, child) in enumerate(geometry)), Val(N))
end

size_scale(geometry::GeometryVectorLike) = isempty(geometry) ? Inf : minimum(size_scale, geometry)
size_scale(geometry::GeometryTuple) = isempty(geometry) ? Inf : minimum(size_scale, geometry)

_geometry_mesh(geometry::GeometryVectorLike; kwargs...) =
    reduce(vcat, (_geometry_mesh(child; kwargs...) for child in geometry); init=Any[])
_geometry_mesh(geometry::GeometryTuple; kwargs...) =
    reduce(vcat, (_geometry_mesh(child; kwargs...) for child in geometry); init=Any[])

end
