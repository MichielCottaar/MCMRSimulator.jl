"""Transformations between physical geometry coordinate systems."""
module Transformations

import StaticArrays: SMatrix, SVector
import ...InternalBoundingBoxes

"""
    Transformation{N, M}

Abstract transformation from an `N`-dimensional coordinate system to an
`M`-dimensional coordinate system.
"""
abstract type Transformation{N, M} end

"""
    Shift{N}(shift)

Translate positions by `shift` in `N` dimensions. Normals are unchanged.
"""
struct Shift{N} <: Transformation{N, N}
    shift::SVector{N, Float64}
end

Shift(shift::AbstractVector{<:Real}) = Shift(SVector{length(shift), Float64}(shift))

"""
    Rotate{N}(matrix)

Apply the orthogonal `N`-dimensional transformation represented by `matrix`.
"""
struct Rotate{N} <: Transformation{N, N}
    matrix::SMatrix{N, N, Float64}
end

function Rotate(matrix::AbstractMatrix{<:Real})
    size(matrix, 1) == size(matrix, 2) ||
        throw(ArgumentError("a Rotate transformation requires a square matrix"))
    N = size(matrix, 1)
    return Rotate{N}(SMatrix{N, N, Float64}(matrix))
end

"""
    Project{N, M}()

Project positions from `N` dimensions to `M` dimensions onto the supported
coordinate axes. Currently, only `3 -> 2` onto the x-y plane and `3 -> 1`
onto the z-axis are supported.
"""
struct Project{N, M} <: Transformation{N, M}
    function Project{N, M}() where {N, M}
        (N, M) in ((3, 2), (3, 1)) ||
            throw(ArgumentError("unsupported Project transformation; only 3 -> 2 and 3 -> 1 are supported"))
        return new{N, M}()
    end
end

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

forward(transformation::Shift, box::InternalBoundingBoxes.InternalBoundingBox) =
    InternalBoundingBoxes.shift(box, transformation.shift)

backward(transformation::Shift, box::InternalBoundingBoxes.InternalBoundingBox) =
    InternalBoundingBoxes.shift(box, -transformation.shift)

forward(transformation::Rotate, position) = transformation.matrix * position
backward(transformation::Rotate, position) = transformation.matrix' * position
forward_normal(transformation::Rotate, normal) = transformation.matrix * normal
backward_normal(transformation::Rotate, normal) = transformation.matrix' * normal

function _rotate_box(transformation::Rotate, box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N}
    center = transformation.matrix * InternalBoundingBoxes._center(box)
    half_size = abs.(transformation.matrix) * InternalBoundingBoxes._half_size(box)
    if InternalBoundingBoxes._is_centered(box)
        return InternalBoundingBoxes.InternalCenteredBoundingRect{N}(half_size)
    end
    return InternalBoundingBoxes.InternalBoundingRect{N}(center, half_size)
end

forward(transformation::Rotate, box::InternalBoundingBoxes.InternalBoundingBox) =
    _rotate_box(transformation, box)

backward(transformation::Rotate, box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N} =
    _rotate_box(Rotate(transformation.matrix'), box)

forward(::Project{3, 2}, position) = SVector(position[1], position[2])
forward(::Project{3, 1}, position) = SVector(position[3])

function forward(::Project{3, 2}, box::InternalBoundingBoxes.InternalBoundingBox)
    half_size = InternalBoundingBoxes._half_size(box)
    if InternalBoundingBoxes._is_cube(box)
        projected_half_size = half_size[1]
    else
        projected_half_size = SVector(half_size[1], half_size[2])
    end
    if InternalBoundingBoxes._is_centered(box)
        if InternalBoundingBoxes._is_cube(box)
            return InternalBoundingBoxes.InternalCenteredBoundingCube{2}(projected_half_size)
        end
        return InternalBoundingBoxes.InternalCenteredBoundingRect{2}(projected_half_size)
    end
    center = InternalBoundingBoxes._center(box)
    if InternalBoundingBoxes._is_cube(box)
        return InternalBoundingBoxes.InternalBoundingCube{2}(SVector(center[1], center[2]), projected_half_size)
    end
    return InternalBoundingBoxes.InternalBoundingRect{2}(SVector(center[1], center[2]), projected_half_size)
end

function forward(::Project{3, 1}, box::InternalBoundingBoxes.InternalBoundingBox)
    half_size = InternalBoundingBoxes._half_size(box)
    projected_half_size = InternalBoundingBoxes._is_cube(box) ? half_size[1] : SVector(half_size[3])
    if InternalBoundingBoxes._is_centered(box)
        if InternalBoundingBoxes._is_cube(box)
            return InternalBoundingBoxes.InternalCenteredBoundingCube{1}(projected_half_size)
        end
        return InternalBoundingBoxes.InternalCenteredBoundingRect{1}(projected_half_size)
    end
    center = InternalBoundingBoxes._center(box)
    if InternalBoundingBoxes._is_cube(box)
        return InternalBoundingBoxes.InternalBoundingCube{1}(SVector(center[3]), projected_half_size)
    end
    return InternalBoundingBoxes.InternalBoundingRect{1}(SVector(center[3]), projected_half_size)
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

end
