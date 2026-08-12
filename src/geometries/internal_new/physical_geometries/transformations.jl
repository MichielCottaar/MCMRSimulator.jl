"""Transformations between physical geometry coordinate systems."""
module Transformations

import StaticArrays: SMatrix, SVector

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

"""Transform a position from the source to the destination coordinate system."""
function forward_position end

"""Transform a position from the destination to the source coordinate system."""
function backward_position end

"""Transform a surface normal from the source to the destination coordinate system."""
function forward_normal end

"""Transform a surface normal from the destination to the source coordinate system."""
function backward_normal end

forward_position(transformation::Shift, position) = position + transformation.shift
backward_position(transformation::Shift, position) = position - transformation.shift
forward_normal(::Shift, normal) = normal
backward_normal(::Shift, normal) = normal

forward_position(transformation::Rotate, position) = transformation.matrix * position
backward_position(transformation::Rotate, position) = transformation.matrix' * position
forward_normal(transformation::Rotate, normal) = transformation.matrix * normal
backward_normal(transformation::Rotate, normal) = transformation.matrix' * normal

forward_position(::Project{3, 2}, position) = SVector(position[1], position[2])
forward_position(::Project{3, 1}, position) = SVector(position[3])

function backward_position(::Project, position)
    throw(ArgumentError("backward_position is not defined for Project transformations"))
end

function forward_normal(::Project, normal)
    throw(ArgumentError("forward_normal is not defined for Project transformations"))
end

function backward_normal(::Project, normal)
    throw(ArgumentError("backward_normal is not defined for Project transformations"))
end

end
