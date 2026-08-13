"""Physical geometries that consist of groups of other physical geometries."""
module Groups

import StaticArrays: SVector
import ...Indices: ObstructionIndex
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside

abstract type GroupGeometry{N} <: PhysicalGeometry{N} end

_is_fully_concrete(::Type{T}) where {T} =
    isconcretetype(T) && all(parameter -> !(parameter isa Type) || _is_fully_concrete(parameter), T.parameters)

struct GeometryVector{N, P<:PhysicalGeometry{N}} <: GroupGeometry{N}
    geometries::Vector{P}

    function GeometryVector{N, P}(geometries::Vector{P}) where {N, P<:PhysicalGeometry{N}}
        _is_fully_concrete(P) || throw(ArgumentError("GeometryVector element types must be fully concrete"))
        new{N, P}(geometries)
    end
end

struct GeometryTuple{N, P<:Tuple{Vararg{PhysicalGeometry{N}}}} <: GroupGeometry{N}
    geometries::P
end

has_inside(::Type{<:GeometryVector{N, P}}) where {N, P} = has_inside(P)
has_inside(::Type{<:GeometryTuple{N, P}}) where {N, P} = any(has_inside, P.parameters)

GeometryVector{N}(geometries::Vector{P}) where {N, P<:PhysicalGeometry{N}} =
    GeometryVector{N, P}(geometries)

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
        child_previous_hit = _previous_hit_for_child(geometry, previous_hit, child_index)
        intersection = detect_intersection(child, start, destination, child_previous_hit)
        if intersection.distance < closest.distance
            closest = _prepend_index(intersection, child_index)
        end
    end
    return closest
end

Base.size(geometry::GeometryVector) = size(geometry.geometries)
Base.axes(geometry::GeometryVector) = axes(geometry.geometries)
Base.length(geometry::GeometryVector) = length(geometry.geometries)
Base.firstindex(geometry::GeometryVector) = firstindex(geometry.geometries)
Base.lastindex(geometry::GeometryVector) = lastindex(geometry.geometries)
Base.getindex(geometry::GeometryVector, index...) = getindex(geometry.geometries, index...)
Base.setindex!(geometry::GeometryVector, value, index...) = setindex!(geometry.geometries, value, index...)
Base.iterate(geometry::GeometryVector, state...) = iterate(geometry.geometries, state...)
Base.eltype(::Type{GeometryVector{N, P}}) where {N, P} = P

Base.length(geometry::GeometryTuple) = length(geometry.geometries)
Base.firstindex(geometry::GeometryTuple) = firstindex(geometry.geometries)
Base.lastindex(geometry::GeometryTuple) = lastindex(geometry.geometries)
Base.getindex(geometry::GeometryTuple, index...) = getindex(geometry.geometries, index...)
Base.iterate(geometry::GeometryTuple, state...) = iterate(geometry.geometries, state...)
Base.eltype(::Type{GeometryTuple{N, P}}) where {N, P} = eltype(P)
Base.Tuple(geometry::GeometryTuple) = geometry.geometries

end
