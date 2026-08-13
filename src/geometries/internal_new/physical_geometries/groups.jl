"""Physical geometries that consist of groups of other physical geometries."""
module Groups

import StaticArrays: SVector
import ...Indices: ObstructionIndex
import ...InternalBoundingBoxes
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, inside_indices, InternalBoundingBox

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

struct GeometryTuple{N, P<:Tuple{Vararg{PhysicalGeometry{N}}}} <: GroupGeometry{N}
    geometries::P
end

has_inside(::Type{<:GeometryVector{N, P}}) where {N, P} = has_inside(P)
has_inside(::Type{<:GeometryVectorBoundingBox{N, P}}) where {N, P} = has_inside(P)
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

function _prepend_index(obstruction_index::ObstructionIndex, index::Int)
    ObstructionIndex(SVector(index, obstruction_index.indices...))
end

function inside_indices(geometry::GeometryVectorLike{N}, position::SVector{N, Float64}) where {N}
    has_inside(typeof(geometry)) || return ObstructionIndex[]
    indices = ObstructionIndex[]
    for (child_index, child) in enumerate(geometry)
        geometry isa GeometryVectorBoundingBox &&
            !InternalBoundingBoxes.isinside(geometry.bounding_boxes[child_index], position) && continue
        for obstruction_index in inside_indices(child, position)
            push!(indices, _prepend_index(obstruction_index, child_index))
        end
    end
    indices
end

function inside_indices(geometry::GeometryTuple{N}, position::SVector{N, Float64}) where {N}
    indices = ObstructionIndex[]
    for (child_index, child) in enumerate(geometry)
        for obstruction_index in inside_indices(child, position)
            push!(indices, _prepend_index(obstruction_index, child_index))
        end
    end
    indices
end

GeometryVector{N}(geometries::Vector{P}) where {N, P<:PhysicalGeometry{N}} =
    GeometryVector{N, P}(geometries)

function GeometryVector(geometries::Vector{P}; bounding_box::Bool=false) where {N, P<:PhysicalGeometry{N}}
    bounding_box ? GeometryVectorBoundingBox(geometries) : GeometryVector{N, P}(geometries)
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

function _prepend_index(intersection::Intersection{N}, index::Int) where {N}
    obstruction_index = intersection.obstruction_index
    indices = SVector(index, obstruction_index.indices...)
    return Intersection(
        intersection.distance,
        intersection.normal,
        intersection.inside,
        ObstructionIndex(indices),
        intersection.hit_gap,
    )
end

function _previous_hit_for_child(
    geometry::GroupGeometry,
    previous_hit::Intersection{3},
    child_index::Int,
)
    Base.isempty(previous_hit) && return Intersection{3}()
    indices = previous_hit.obstruction_index.indices
    isempty(indices) && throw(ArgumentError("a non-empty previous hit must have an obstruction index"))
    indices[1] == child_index || return Intersection{3}()
    remaining_indices = SVector{length(indices) - 1, Int}(indices[2:end])
    return Intersection(
        previous_hit.distance,
        previous_hit.normal,
        previous_hit.inside,
        ObstructionIndex(remaining_indices),
        previous_hit.hit_gap,
    )
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
            closest = _prepend_index(intersection, child_index)
        end
    end
    return closest
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

end
