"""Compact internal representations of axis-aligned bounding boxes."""
module InternalBoundingBoxes

import StaticArrays: SVector

abstract type InternalBoundingBox{N} end

struct InternalBoundingCube{N} <: InternalBoundingBox{N}
    center::SVector{N, Float64}
    half_size::Float64
end

struct InternalCenteredBoundingCube{N} <: InternalBoundingBox{N}
    half_size::Float64
end

struct InternalBoundingRect{N} <: InternalBoundingBox{N}
    center::SVector{N, Float64}
    half_size::SVector{N, Float64}
end

struct InternalCenteredBoundingRect{N} <: InternalBoundingBox{N}
    half_size::SVector{N, Float64}
end

_is_centered(::InternalCenteredBoundingCube) = true
_is_centered(::InternalBoundingCube) = false
_is_centered(::InternalCenteredBoundingRect) = true
_is_centered(::InternalBoundingRect) = false

_is_cube(::InternalCenteredBoundingCube) = true
_is_cube(::InternalBoundingCube) = true
_is_cube(::InternalCenteredBoundingRect) = false
_is_cube(::InternalBoundingRect) = false

InternalBoundingCube(center::AbstractVector{<:Real}, half_size::Real) =
    InternalBoundingCube{length(center)}(SVector{length(center), Float64}(center), Float64(half_size))

InternalBoundingRect(center::AbstractVector{<:Real}, half_size::AbstractVector{<:Real}) =
    InternalBoundingRect{length(center)}(
        SVector{length(center), Float64}(center),
        SVector{length(center), Float64}(half_size),
    )

shift(box::InternalCenteredBoundingCube{N}, displacement::SVector{N, Float64}) where {N} =
    InternalBoundingCube{N}(displacement, box.half_size)

shift(box::InternalBoundingCube{N}, displacement::SVector{N, Float64}) where {N} =
    InternalBoundingCube{N}(box.center + displacement, box.half_size)

shift(box::InternalCenteredBoundingRect{N}, displacement::SVector{N, Float64}) where {N} =
    InternalBoundingRect{N}(displacement, box.half_size)

shift(box::InternalBoundingRect{N}, displacement::SVector{N, Float64}) where {N} =
    InternalBoundingRect{N}(box.center + displacement, box.half_size)

shift(box::InternalBoundingBox{N}, displacement::AbstractVector{<:Real}) where {N} =
    shift(box, SVector{N, Float64}(displacement))

function InternalBoundingBox(half_size::Real, center=nothing)
    if isnothing(center)
        return InternalCenteredBoundingCube{3}(Float64(half_size))
    end
    return InternalBoundingCube(center, half_size)
end

function InternalBoundingBox(half_size::AbstractVector{<:Real}, center=nothing)
    if isnothing(center)
        return InternalCenteredBoundingRect{length(half_size)}(
            SVector{length(half_size), Float64}(half_size),
        )
    end
    return InternalBoundingRect(center, half_size)
end

lower(box::InternalBoundingBox) = _center(box) - _half_size(box)
upper(box::InternalBoundingBox) = _center(box) + _half_size(box)

_center(box::InternalBoundingBox{N}) where {N} = _is_centered(box) ?
    zero(SVector{N, Float64}) : _stored_center(box)

_stored_center(box::InternalBoundingCube) = box.center
_stored_center(box::InternalBoundingRect) = box.center

_half_size(box::InternalBoundingBox{N}) where {N} = _is_cube(box) ?
    SVector{N, Float64}(fill(box.half_size, N)) : box.half_size

isinside(box::InternalBoundingBox, position::AbstractVector) =
    all(position .>= lower(box)) && all(position .<= upper(box))

function could_intersect(box::InternalBoundingBox{N}, start::AbstractVector, destination::AbstractVector) where {N}
    lower_bound = lower(box)
    upper_bound = upper(box)
    return all(
        (start[i] <= upper_bound[i] || destination[i] <= upper_bound[i]) &&
        (start[i] >= lower_bound[i] || destination[i] >= lower_bound[i])
        for i in 1:N
    )
end

function does_intersect(box::InternalBoundingBox{N}, start::AbstractVector, destination::AbstractVector) where {N}
    lower_bound = lower(box)
    upper_bound = upper(box)
    t_min = 0.0
    t_max = 1.0

    for i in 1:N
        displacement = destination[i] - start[i]
        if iszero(displacement)
            isinside_interval = lower_bound[i] <= start[i] <= upper_bound[i]
            isinside_interval || return false
        else
            t1 = (lower_bound[i] - start[i]) / displacement
            t2 = (upper_bound[i] - start[i]) / displacement
            t_min = max(t_min, min(t1, t2))
            t_max = min(t_max, max(t1, t2))
            t_min <= t_max || return false
        end
    end
    return true
end

function corners(box::InternalBoundingBox{N}) where {N}
    lower_bound = lower(box)
    upper_bound = upper(box)
    return [
        SVector{N, Float64}(ntuple(i -> (mask >> (i - 1)) & 1 == 0 ? lower_bound[i] : upper_bound[i], N))
        for mask in 0:(2^N - 1)
    ]
end

end
