"""Transformations between physical geometry coordinate systems."""
module Transformations

import StaticArrays: SMatrix, SVector
import LinearAlgebra: norm
import ...InternalBoundingBoxes
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, inside_indices, InternalBoundingBox

"""
    Transformation{N, M, P} <: PhysicalGeometry{N}

Transformation of a geometry from an `N`-dimensional coordinate system to an
`M`-dimensional coordinate system.

`forward` transforms map from the global N-dimensional space to the local geometry-specific M-dimensional space.
`backward` are the inverse transformations.
"""
abstract type Transformation{N, M, P<:PhysicalGeometry{M}} <: PhysicalGeometry{N} end

has_inside(::Type{<:Transformation{N, M, P}}) where {N, M, P} = has_inside(P)

InternalBoundingBox(transformation::Transformation) =
    backward(transformation, InternalBoundingBox(transformation.geometry))

inside_indices(transformation::Transformation, position) =
    inside_indices(transformation.geometry, forward(transformation, position))

"""
    Shift{N}(shift)

Translate positions by `shift` in `N` dimensions. Normals are unchanged.
"""
struct Shift{N, P<:PhysicalGeometry{N}} <: Transformation{N, N, P}
    geometry::P
    shift::SVector{N, Float64}
end

Shift(geometry::P, shift::AbstractVector{<:Real}) where {N, P<:PhysicalGeometry{N}} =
    Shift{N, P}(geometry, SVector{N, Float64}(shift))

"""
    Scale{N}(scale)

Uniformly scale positions by a positive `scale` in `N` dimensions.
"""
struct Scale{N, P<:PhysicalGeometry{N}} <: Transformation{N, N, P}
    geometry::P
    scale::Float64

    function Scale{N, P}(geometry::P, scale::Real) where {N, P<:PhysicalGeometry{N}}
        scale > 0 || throw(ArgumentError("a Scale transformation requires a positive scale"))
        scale = Float64(scale)
        return new{N, P}(geometry, scale)
    end
end

Scale(geometry::P, scale::Real) where {N, P<:PhysicalGeometry{N}} = Scale{N, P}(geometry, scale)

"""
    Rotate(geometry, matrix)

Map coordinates from the outer `N`-dimensional space to the child `M`-dimensional
space using an orthonormal matrix. For `M == N` this is a rotation; for `M < N`
it is a rotation followed by projection.
"""
struct Rotate{N, M, P<:PhysicalGeometry{M}} <: Transformation{N, M, P}
    geometry::P
    matrix::SMatrix{M, N, Float64}

    function Rotate{N, M, P}(geometry::P, matrix::SMatrix{M, N, Float64}) where {N, M, P<:PhysicalGeometry{M}}
        _validate_rotate_matrix(matrix, N, M)
        new{N, M, P}(geometry, matrix)
    end
end

function _validate_rotate_matrix(matrix, N, M)
    M <= N || throw(ArgumentError("a Rotate transformation cannot increase dimension"))
    all(
        isapprox(
            sum(matrix[i, k] * matrix[j, k] for k in 1:N),
            i == j ? 1.0 : 0.0;
            atol=1e-10,
            rtol=1e-10,
        )
        for i in 1:M, j in 1:M
    ) || throw(ArgumentError("Rotate matrix rows must be orthonormal"))
    nothing
end

function Rotate(geometry::P, matrix::AbstractMatrix{<:Real}) where {M, P<:PhysicalGeometry{M}}
    N = size(matrix, 2)
    size(matrix, 1) == M || throw(ArgumentError("rotation matrix rows must match geometry dimension"))
    matrix = SMatrix{M, N, Float64}(matrix)
    Rotate{N, M, P}(geometry, matrix)
end

InternalBoundingBox(transformation::Rotate{N, M}) where {N, M} =
    N == M ? backward(transformation, InternalBoundingBox(transformation.geometry)) :
    throw(ArgumentError("dimension-reducing Rotate transformations do not have a finite bounding box"))

"""Transform a value from the source to the destination coordinate system."""
function forward end

"""Transform a value from the destination to the source coordinate system."""
function backward end

"""Transform a surface normal from the source to the destination coordinate system."""
function forward_normal end

"""Transform a surface normal from the destination to the source coordinate system."""
function backward_normal end

forward(transformation::Shift, position) = position + transformation.shift
backward(transformation::Shift, position) = position - transformation.shift
forward_normal(::Shift, normal) = normal
backward_normal(::Shift, normal) = normal

forward(transformation::Scale, position) = transformation.scale * position
backward(transformation::Scale, position) = position / transformation.scale
forward_normal(::Scale, normal) = normal
backward_normal(::Scale, normal) = normal

forward(transformation::Shift, box::InternalBoundingBoxes.InternalBoundingBox) =
    InternalBoundingBoxes.shift(box, transformation.shift)

backward(transformation::Shift, box::InternalBoundingBoxes.InternalBoundingBox) =
    InternalBoundingBoxes.shift(box, -transformation.shift)

function _box_center(box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N}
    return InternalBoundingBoxes.center(box)
end

function _box_half_size(box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N}
    return InternalBoundingBoxes.half_size(box)
end

function _scale_box(transformation::Scale, box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N}
    half_size = transformation.scale .* _box_half_size(box)
    center = transformation.scale .* _box_center(box)
    return InternalBoundingBoxes.InternalBoundingBox{N}(
        half_size,
        center,
    )
end

forward(transformation::Scale, box::InternalBoundingBoxes.InternalBoundingBox) =
    _scale_box(transformation, box)

backward(transformation::Scale, box::InternalBoundingBoxes.InternalBoundingBox) =
    _scale_box(Scale(transformation.geometry, inv(transformation.scale)), box)

forward(transformation::Rotate, position) = transformation.matrix * position
backward(transformation::Rotate{N, M}, position) where {N, M} =
    N == M ? transformation.matrix' * position : throw(ArgumentError("backward is not defined for dimension-reducing Rotate transformations"))
forward_normal(transformation::Rotate{N, M}, normal) where {N, M} =
    N == M ? transformation.matrix * normal : throw(ArgumentError("forward_normal is not defined for dimension-reducing Rotate transformations"))
backward_normal(transformation::Rotate, normal) = begin
    result = transformation.matrix' * normal
    result ./ norm(result)
end

function _rotate_box(transformation::Rotate{N, M}, box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N, M}
    center = transformation.matrix * _box_center(box)
    half_size = abs.(transformation.matrix) * _box_half_size(box)
    return InternalBoundingBoxes.InternalBoundingBox{M}(half_size, center)
end

forward(transformation::Rotate, box::InternalBoundingBoxes.InternalBoundingBox) =
    _rotate_box(transformation, box)

backward(transformation::Rotate{N, M}, box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N, M} =
    N == M ? InternalBoundingBoxes.InternalBoundingBox{N}(
        abs.(transformation.matrix') * _box_half_size(box),
        transformation.matrix' * _box_center(box),
    ) : throw(ArgumentError("backward is not defined for dimension-reducing Rotate transformations"))

function detect_intersection(
    transformation::Transformation{N, M, P},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
) where {N, M, P}
    child_intersection = detect_intersection(
        transformation.geometry,
        forward(transformation, start),
        forward(transformation, destination),
        previous_hit,
    )
    Base.isempty(child_intersection) && return Intersection{N}()
    return Intersection(
        child_intersection.distance,
        backward_normal(transformation, child_intersection.normal),
        child_intersection.inside,
        child_intersection.obstruction_index,
        child_intersection.hit_gap,
    )
end

end
