"""Transformations between physical geometry coordinate systems."""
module Transformations

import StaticArrays: SMatrix, SVector
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
    Rotate{N}(matrix)

Apply the orthogonal `N`-dimensional transformation represented by `matrix`.
"""
struct Rotate{N, P<:PhysicalGeometry{N}} <: Transformation{N, N, P}
    geometry::P
    matrix::SMatrix{N, N, Float64}
end

function Rotate(geometry::P, matrix::AbstractMatrix{<:Real}) where {N, P<:PhysicalGeometry{N}}
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("a Rotate transformation requires a square matrix"))
    size(matrix, 1) == N || throw(ArgumentError("rotation dimension must match geometry dimension"))
    return Rotate{N, P}(geometry, SMatrix{N, N, Float64}(matrix))
end

"""
    Project{N, M}()

Project positions from `N` dimensions to `M` dimensions onto the supported
coordinate axes. Currently, only `3 -> 2` onto the x-y plane and `3 -> 1`
onto the z-axis are supported.
"""
struct Project{N, M, P<:PhysicalGeometry{M}} <: Transformation{N, M, P}
    geometry::P

    function Project{N, M, P}(geometry::P) where {N, M, P<:PhysicalGeometry{M}}
        (N, M) in ((3, 2), (3, 1)) ||
            throw(ArgumentError("unsupported Project transformation; only 3 -> 2 and 3 -> 1 are supported"))
        return new{N, M, P}(geometry)
    end
end

Project{N, M}(geometry::P) where {N, M, P<:PhysicalGeometry{M}} = Project{N, M, P}(geometry)

InternalBoundingBox(::Project) = throw(ArgumentError("Project transformations do not have a finite backward bounding box"))

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
backward(transformation::Rotate, position) = transformation.matrix' * position
forward_normal(transformation::Rotate, normal) = transformation.matrix * normal
backward_normal(transformation::Rotate, normal) = transformation.matrix' * normal

function _rotate_box(transformation::Rotate, box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N}
    center = transformation.matrix * _box_center(box)
    half_size = abs.(transformation.matrix) * _box_half_size(box)
    return InternalBoundingBoxes.InternalBoundingBox{N}(half_size, center)
end

forward(transformation::Rotate, box::InternalBoundingBoxes.InternalBoundingBox) =
    _rotate_box(transformation, box)

backward(transformation::Rotate, box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N} =
    _rotate_box(Rotate(transformation.geometry, transformation.matrix'), box)

forward(::Project{3, 2}, position) = SVector(position[1], position[2])
forward(::Project{3, 1}, position) = SVector(position[3])

function forward(::Project{3, 2}, box::InternalBoundingBoxes.InternalBoundingBox)
    half_size = _box_half_size(box)
    projected_half_size = SVector(half_size[1], half_size[2])
    center = _box_center(box)
    projected_center = SVector(center[1], center[2])
    return InternalBoundingBoxes.InternalBoundingBox{2}(projected_half_size, projected_center)
end

function forward(::Project{3, 1}, box::InternalBoundingBoxes.InternalBoundingBox)
    half_size = _box_half_size(box)
    projected_half_size = SVector(half_size[3])
    center = _box_center(box)
    projected_center = SVector(center[3])
    return InternalBoundingBoxes.InternalBoundingBox{1}(projected_half_size, projected_center)
end

function backward(::Project, position)
    throw(ArgumentError("backward is not defined for Project transformations"))
end

function backward(::Project, ::InternalBoundingBoxes.InternalBoundingBox)
    throw(ArgumentError("backward is not defined for Project transformations"))
end

function forward_normal(::Project, normal)
    throw(ArgumentError("forward_normal is not defined for Project transformations"))
end

function backward_normal(::Project, normal)
    throw(ArgumentError("backward_normal is not defined for Project transformations"))
end

backward_normal(::Project{3, 2}, normal) = SVector(normal[1], normal[2], 0.0)
backward_normal(::Project{3, 1}, normal) = SVector(0.0, 0.0, normal[1])

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
