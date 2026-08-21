"""Physical geometries that consist of groups of other physical geometries."""
module Groups

import StaticArrays: SVector
import ...InternalBoundingBoxes
import ..GridDispatch: IntersectionGrid, GridIterator
import ..PhysicalGeometries: PhysicalGeometry, child_type, find_intersection, get_child, has_inside, has_single_inside, inside_indices_eltype, isinside_single, inside_indices, InternalBoundingBox
import ..PhysicalGeometries: random_surface_positions, size_scale, distance_to_surface, _geometry_mesh
import ...Properties: GeometryProperties, GeometryLeafProperties, GeometryVectorProperties, GeometryTupleProperties

abstract type GroupGeometry{N, P} <: PhysicalGeometry{N} end
abstract type GeometryVectorLike{N, P<:PhysicalGeometry{N}} <: GroupGeometry{N, P} end

function Base.show(io::IO, ::Type{T}) where {N, P, T <: GeometryVectorLike{N, P}}
    print(io, "GeometryVector{")
    show(io, P)
    print(io, "}")
end

group_geometries(group::GeometryVectorLike; include_gap=true) = group.geometries

function distance_to_surface(group::GroupGeometry{N}, position::SVector{N, Float64}) where {N}
    minimum(
        (distance_to_surface(child, position)
         for child in group_geometries(group; include_gap=false));
        init=Inf,
    )
end

child_type(::Type{<:GroupGeometry{N, P}}) where {N, P} = P

function _prepend_type(::Type{Prefix}, ::Type{T}) where {Prefix, T}
    T === Union{} && return Union{}
    T isa Union && return Union{
        (_prepend_type(Prefix, element) for element in Base.uniontypes(T))...,
    }
    T <: Tuple && return Tuple{Prefix, T.parameters...}
    throw(MethodError(_prepend_type, (Type{Prefix}, Type{T})))
end

function inside_indices_eltype(::Type{T}) where {T}
    T isa Union || throw(MethodError(inside_indices_eltype, (Type{T},)))
    Union{(inside_indices_eltype(element) for element in Base.uniontypes(T))...}
end

inside_indices_eltype(::Type{<:GroupGeometry{N, P}}) where {N, P} =
    _prepend_type(Int, inside_indices_eltype(P))


function find_intersection(group::GroupGeometry{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, previous_hit=nothing) where {N}
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

group_geometries(group::GeometryTuple; include_gap=true) = group.geometries

function get_child(group::GeometryVectorLike, indices::Tuple)
    child_index = indices[1]
    group_geometries(group)[child_index], indices[2:end]
end

function get_child(group::GeometryTuple, indices::Tuple)
    child_index = indices[1]
    group_geometries(group)[child_index], indices[2:end]
end

child_type(::Type{<:GeometryTuple{N, P}}) where {N, P} = Union{P.parameters...}

inside_indices_eltype(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    _prepend_type(Int, inside_indices_eltype(Union{P.parameters...}))

function inside_candidates end

inside_candidates(group::GeometryVector, position) = enumerate(group_geometries(group))

inside_candidates(group::GeometryVectorBoundingBox, position) =
    (
        (index, group_geometries(group)[index])
        for index in eachindex(group_geometries(group))
        if InternalBoundingBoxes.isinside(group.bounding_boxes[index], position)
    )

function inside_candidates(group::GeometryVectorGrid{N}, position::SVector{N, Float64}) where {N}
    coordinate = _grid_coordinate(group, position)
    isnothing(coordinate) && return ()
    ((index, group_geometries(group)[index]) for index in group.grid.indices[coordinate...])
end

inside_candidates(group::GeometryTuple, position) = enumerate(group_geometries(group))

intersection_candidates(group::Union{GeometryVector, GeometryTuple}, start, destination) =
    ((index, child, 0.0) for (index, child) in enumerate(group_geometries(group)))

intersection_candidates(group::GeometryVectorBoundingBox, start, destination) =
    (
        (index, group_geometries(group)[index], 0.0)
        for index in eachindex(group_geometries(group))
        if InternalBoundingBoxes.does_intersect(group.bounding_boxes[index], start, destination)
    )

function intersection_candidates(group::GeometryVectorGrid, start, destination)
    start_coordinate = _grid_coordinate(group, start)
    destination_coordinate = _grid_coordinate(group, destination)
    if !isnothing(start_coordinate) && start_coordinate == destination_coordinate
        return (
            (index, group_geometries(group)[index], 0.0)
            for index in group.grid.indices[start_coordinate...]
        )
    end
    (
        (child_index, group_geometries(group)[child_index], dist_all_checked)
        for (child_index, dist_all_checked) in GridIterator(group.grid, start, destination)
    )
end

has_inside(::Type{<:GeometryVectorLike{N, P}}) where {N, P} = has_inside(P)
has_inside(::Type{<:GeometryTuple{N, P}}) where {N, P} = any(has_inside, P.parameters)
has_single_inside(::Type{<:GroupGeometry}) = false

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

function inside_indices_for_any_type(geometry, position, intersection)
    geometry_type = typeof(geometry)
    has_inside(geometry_type) || return Tuple{}[]
    if has_single_inside(geometry_type)
        return isinside_single(geometry, position, intersection) ?
            [()] : Tuple{}[]
    end
    return inside_indices(geometry, position, intersection)
end

function inside_indices(
    geometry::GroupGeometry{N, P},
    position::SVector{N, Float64},
    intersection=nothing,
) where {N, P}
    indices = Vector{inside_indices_eltype(typeof(geometry))}()
    if !has_inside(P)
        return indices
    end
    single_inside = has_single_inside(P)
    for (child_index, child) in inside_candidates(geometry, position)
        child_intersection = !isnothing(intersection) && intersection[1] == child_index ? intersection[2:end] : nothing
        if single_inside
            if isinside_single(child, position, child_intersection)
                push!(indices, (child_index, ))
            end
        else
            append!(indices, [(child_index, new_indices...) for new_indices in inside_indices(child, position, child_intersection)])
        end
    end
    return indices
end

function _grid_coordinate(
    geometry::GeometryVectorGrid{N},
    position::SVector{N, Float64},
) where {N}
    coordinate = Int.(floor.((position - InternalBoundingBoxes.lower(geometry.grid.grid_bounding_box)) .* geometry.grid.inv_resolution)) .+ 1
    any(coordinate .< 1) || any(coordinate .> size(geometry.grid.indices)) ? nothing : SVector{N, Int}(coordinate)
end

function inside_indices(
    geometry::GeometryTuple{N},
    position::SVector{N, Float64},
    intersection=nothing
) where {N}
    indices = Vector{inside_indices_eltype(typeof(geometry))}()
    for (child_index, child) in enumerate(geometry)
        child_intersection = !isnothing(intersection) && intersection[1] == child_index ? intersection[2:end] : nothing
        append!(indices, [(child_index, new_indices...) for new_indices in inside_indices_for_any_type(child, position, child_intersection)])
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

function Base.show(io::IO, ::Type{T}) where {N, P, T <: GeometryTuple{N, P}}
    print(io, "GeometryTuple{")
    for (index, child) in enumerate(P.parameters)
        index > 1 && print(io, ", ")
        show(io, child)
    end
    print(io, "}")
end

_density_child(density::GeometryLeafProperties, ::Int) = density
_density_child(density::GeometryVectorProperties, index::Int) = density.properties[index]
_density_child(density::GeometryTupleProperties, index::Int) = density.properties[index]
_density_child(density::GeometryProperties, ::SVector) = density

function _combine(draws, ::Val{N}) where {N}
    draws = collect(draws)
    positions = reduce(vcat, (draw[1] for draw in draws); init=SVector{N, Float64}[])
    indices = reduce(vcat, (draw[2] for draw in draws); init=Tuple[])
    positions, indices
end

function random_surface_positions(geometry::GroupGeometry{N}, density::GeometryProperties,
    bounding_box::InternalBoundingBox{N}, scale_density) where {N}
    _combine((let
        values = random_surface_positions(child, _density_child(density, index), bounding_box, scale_density)
        (values[1], [(index, child_index...) for child_index in values[2]])
    end for (index, child) in enumerate(group_geometries(geometry; include_gap=false))), Val(N))
end

size_scale(geometry::GeometryVectorLike) = isempty(geometry) ? Inf : minimum(size_scale, geometry)
size_scale(geometry::GeometryTuple) = isempty(geometry) ? Inf : minimum(size_scale, geometry)

_geometry_mesh(geometry::GeometryVectorLike; kwargs...) =
    reduce(vcat, (_geometry_mesh(child; kwargs...) for child in geometry); init=Any[])
_geometry_mesh(geometry::GeometryTuple; kwargs...) =
    reduce(vcat, (_geometry_mesh(child; kwargs...) for child in geometry); init=Any[])

end
