"""Physical geometries that consist of groups of other physical geometries."""
module Groups

import StaticArrays: SVector
import ...Indices: ObstructionIndex, add_index, remove_index
import ...InternalBoundingBoxes
import ..GridDispatch: detect_intersection_grid
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, inside_indices, InternalBoundingBox
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh
import ...Properties: GeometryProperties, GeometryLeafProperties, GeometryVectorProperties, GeometryTupleProperties

abstract type GroupGeometry{N} <: PhysicalGeometry{N} end
abstract type GeometryVectorLike{N, P<:PhysicalGeometry{N}} <: GroupGeometry{N} end

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
    bounding_box::InternalBoundingBoxes.InternalBoundingBox{N}
    inv_resolution::SVector{N, Float64}
    indices::Array{Vector{Int}, N}

    function GeometryVectorGrid{N, P}(
        geometries::Vector{P},
        bounding_box::InternalBoundingBoxes.InternalBoundingBox{N},
        inv_resolution::SVector{N, Float64},
        indices::Array{Vector{Int}, N},
    ) where {N, P<:PhysicalGeometry{N}}
        _is_fully_concrete(P) || throw(ArgumentError("GeometryVector element types must be fully concrete"))
        new{N, P}(geometries, bounding_box, inv_resolution, indices)
    end
end

struct GeometryTuple{N, P<:Tuple{Vararg{PhysicalGeometry{N}}}} <: GroupGeometry{N}
    geometries::P
end

has_inside(::Type{<:GeometryVector{N, P}}) where {N, P} = has_inside(P)
has_inside(::Type{<:GeometryVectorBoundingBox{N, P}}) where {N, P} = has_inside(P)
has_inside(::Type{<:GeometryVectorGrid{N, P}}) where {N, P} = has_inside(P)
has_inside(::Type{<:GeometryTuple{N, P}}) where {N, P} = any(has_inside, P.parameters)

function InternalBoundingBox(geometry::GroupGeometry{N}) where {N}
    isempty(geometry) && throw(ArgumentError("cannot construct a bounding box for an empty geometry group"))
    boxes = [InternalBoundingBoxes.InternalBoundingBox(child) for child in geometry]
    lower_bound = reduce((a, b) -> min.(a, b), (InternalBoundingBoxes.lower(box) for box in boxes))
    upper_bound = reduce((a, b) -> max.(a, b), (InternalBoundingBoxes.upper(box) for box in boxes))
    InternalBoundingBoxes.InternalBoundingBox(
        (upper_bound - lower_bound) / 2,
        (upper_bound + lower_bound) / 2,
    )
end

InternalBoundingBox(geometry::GeometryVectorGrid) = geometry.bounding_box

function _inside_child_intersection(
    intersection::Intersection{N},
    child_index::Int,
) where {N}
    Base.isempty(intersection) && return Intersection{N}()
    isempty(intersection.obstruction_index.indices) && return Intersection{N}()
    index, child_intersection = remove_index(intersection)
    index == child_index ? child_intersection : Intersection{N}()
end

function inside_indices(
    geometry::GeometryVectorLike{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N}
    has_inside(typeof(geometry)) || return ObstructionIndex[]
    indices = ObstructionIndex[]
    for (child_index, child) in enumerate(geometry)
        geometry isa GeometryVectorBoundingBox &&
            !InternalBoundingBoxes.isinside(geometry.bounding_boxes[child_index], position) && continue
        child_intersection = _inside_child_intersection(intersection, child_index)
        for obstruction_index in inside_indices(child, position, child_intersection)
            push!(indices, add_index(obstruction_index, child_index))
        end
    end
    indices
end

function _grid_coordinate(
    geometry::GeometryVectorGrid{N},
    position::SVector{N, Float64},
) where {N}
    coordinate = Int.(floor.((position - InternalBoundingBoxes.lower(geometry.bounding_box)) .* geometry.inv_resolution)) .+ 1
    any(coordinate .< 1) || any(coordinate .> size(geometry.indices)) ? nothing : SVector{N, Int}(coordinate)
end

function inside_indices(
    geometry::GeometryVectorGrid{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N}
    has_inside(typeof(geometry)) || return ObstructionIndex[]
    coordinate = _grid_coordinate(geometry, position)
    isnothing(coordinate) && return ObstructionIndex[]
    indices = ObstructionIndex[]
    for child_index in geometry.indices[coordinate...]
        child_intersection = _inside_child_intersection(intersection, child_index)
        for obstruction_index in inside_indices(geometry.geometries[child_index], position, child_intersection)
            push!(indices, add_index(obstruction_index, child_index))
        end
    end
    indices
end

function _grid_candidates(geometry::GeometryVectorGrid, coordinate)
    (any(coordinate .< 1) || any(coordinate .> size(geometry.indices))) && return Int[]
    geometry.indices[coordinate...]
end

function inside_indices(
    geometry::GeometryTuple{N},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N}
    indices = ObstructionIndex[]
    for (child_index, child) in enumerate(geometry)
        child_intersection = _inside_child_intersection(intersection, child_index)
        for obstruction_index in inside_indices(child, position, child_intersection)
            push!(indices, add_index(obstruction_index, child_index))
        end
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

function _grid_dimensions(box::InternalBoundingBoxes.InternalBoundingBox{N}, resolution) where {N}
    extent = InternalBoundingBoxes.upper(box) - InternalBoundingBoxes.lower(box)
    if isinf(resolution)
        return SVector{N, Int}(ones(Int, N)), extent
    end
    resolution > 0 || throw(ArgumentError("grid_resolution must be positive"))
    dimensions = max.(SVector{N, Int}(ceil.(Int, extent ./ resolution)), 1)
    return dimensions, extent ./ dimensions
end

function GeometryVectorGrid(
    geometries::Vector{P};
    grid_resolution=nothing,
) where {N, P<:PhysicalGeometry{N}}
    isempty(geometries) && throw(ArgumentError("cannot construct a grid for an empty geometry vector"))
    boxes = InternalBoundingBox.(geometries)
    lower_bound = reduce((a, b) -> min.(a, b), (InternalBoundingBoxes.lower(box) for box in boxes))
    upper_bound = reduce((a, b) -> max.(a, b), (InternalBoundingBoxes.upper(box) for box in boxes))
    box = InternalBoundingBoxes.InternalBoundingBox(
        (upper_bound - lower_bound) / 2,
        (upper_bound + lower_bound) / 2,
    )
    extent = upper_bound - lower_bound
    resolution = isnothing(grid_resolution) ?
        (maximum(extent) / max(ceil(Int, length(geometries)^(1 / N)), 1)) : grid_resolution
    dimensions, cell_size = _grid_dimensions(box, resolution)
    actual_box = InternalBoundingBoxes.InternalBoundingBox(
        cell_size .* dimensions / 2,
        lower_bound + cell_size .* dimensions / 2,
    )
    inv_resolution = SVector{N, Float64}(1 ./ cell_size)
    indices = InternalBoundingBoxes.grid_indices(actual_box, dimensions, boxes)
    GeometryVectorGrid{N, P}(geometries, actual_box, inv_resolution, indices)
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

function _previous_hit_for_child(
    geometry::GroupGeometry,
    previous_hit::Intersection{3},
    child_index::Int,
)
    Base.isempty(previous_hit) && return Intersection{3}()
    isempty(previous_hit.obstruction_index.indices) &&
        throw(ArgumentError("a non-empty previous hit must have an obstruction index"))
    index, child_previous_hit = remove_index(previous_hit)
    index == child_index ? child_previous_hit : Intersection{3}()
end

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
    geometry::GroupGeometry{N},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
) where {N}
    if !Base.isempty(previous_hit)
        indices = previous_hit.obstruction_index.indices
        isempty(indices) && throw(ArgumentError("a non-empty previous hit must have an obstruction index"))
        1 <= indices[1] <= length(geometry) ||
            throw(ArgumentError("previous-hit obstruction index is not a child of this group"))
    end

    closest = Intersection{N}()
    for (child_index, child) in enumerate(geometry)
        _could_intersect(geometry, child_index, start, destination) || continue
        child_previous_hit = _previous_hit_for_child(geometry, previous_hit, child_index)
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
    previous_hit::Intersection{3}=Intersection{3}(),
) where {N}
    return detect_intersection_grid(
        geometry.bounding_box,
        geometry.inv_resolution,
        geometry.indices,
        start,
        destination,
        previous_hit,
    ) do child_index, previous
        child_previous_hit = _previous_hit_for_child(geometry, previous, child_index)
        intersection = detect_intersection(geometry.geometries[child_index], start, destination, child_previous_hit)
        add_index(intersection, child_index)
    end
end

Base.size(geometry::GeometryVectorLike) = size(geometry.geometries)
Base.axes(geometry::GeometryVectorLike) = axes(geometry.geometries)
Base.length(geometry::GeometryVectorLike) = length(geometry.geometries)
Base.firstindex(geometry::GeometryVectorLike) = firstindex(geometry.geometries)
Base.lastindex(geometry::GeometryVectorLike) = lastindex(geometry.geometries)
Base.getindex(geometry::GeometryVectorLike, index...) = getindex(geometry.geometries, index...)
Base.setindex!(geometry::GeometryVectorLike, value, index...) = setindex!(geometry.geometries, value, index...)
Base.iterate(geometry::GeometryVectorLike, state...) = iterate(geometry.geometries, state...)
Base.eltype(::Type{<:GeometryVectorLike{N, P}}) where {N, P} = P

Base.length(geometry::GeometryTuple) = length(geometry.geometries)
Base.firstindex(geometry::GeometryTuple) = firstindex(geometry.geometries)
Base.lastindex(geometry::GeometryTuple) = lastindex(geometry.geometries)
Base.getindex(geometry::GeometryTuple, index...) = getindex(geometry.geometries, index...)
Base.iterate(geometry::GeometryTuple, state...) = iterate(geometry.geometries, state...)
Base.eltype(::Type{GeometryTuple{N, P}}) where {N, P} = eltype(P)
Base.Tuple(geometry::GeometryTuple) = geometry.geometries

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
