"""Transformations between physical geometry coordinate systems."""
module Transformations

import StaticArrays: SMatrix, SVector
import LinearAlgebra: norm, nullspace, cross
import Random: rand
import ...InternalBoundingBoxes
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, has_single_inside, isinside_single, inside_indices, InternalBoundingBox
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh, _mesh_result, _translate_native
import ...Properties: GeometryProperties

"""
    Transformation{N, M, P} <: PhysicalGeometry{N}

Transformation of a geometry from an `N`-dimensional coordinate system to an
`M`-dimensional coordinate system.

`forward` transforms map from the global N-dimensional space to the local geometry-specific M-dimensional space.
`backward` are the inverse transformations.
"""
abstract type Transformation{N, M, P<:PhysicalGeometry{M}} <: PhysicalGeometry{N} end

has_inside(::Type{<:Transformation{N, M, P}}) where {N, M, P} = has_inside(P)
has_single_inside(::Type{<:Transformation{N, M, P}}) where {N, M, P} = has_single_inside(P)

InternalBoundingBox(transformation::Transformation) =
    backward(transformation, InternalBoundingBox(transformation.geometry))

function isinside_single(
    transformation::Transformation{N, M},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N, M}
    isinside_single(
        transformation.geometry,
        forward(transformation, position),
        intersection,
    )
end

function inside_indices(
    transformation::Transformation{N, M},
    position::SVector{N, Float64},
    intersection::Intersection{N}=Intersection{N}(),
) where {N, M}
    inside_indices(
        transformation.geometry,
        forward(transformation, position),
        intersection,
    )
end

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

size_scale(geometry::Shift) = size_scale(geometry.geometry)
size_scale(geometry::Rotate) = size_scale(geometry.geometry)
size_scale(geometry::Scale) = geometry.scale * size_scale(geometry.geometry)

function _deproject_positions(transformation::Rotate{N, M}, positions,
    bounding_box::InternalBoundingBox{N}) where {N, M}
    basis = nullspace(Matrix(transformation.matrix))
    center = InternalBoundingBoxes.center(bounding_box)
    half_size = InternalBoundingBoxes.half_size(bounding_box)
    null_center = basis' * center
    null_half_size = abs.(basis)' * half_size
    [transformation.matrix' * position + basis * (null_center +
        (2 .* rand(length(null_center)) .- 1) .* null_half_size) for position in positions]
end

_projected_scale(transformation::Rotate{N, M}, bounding_box::InternalBoundingBox{N}) where {N, M} =
    prod(2 .* (abs.(nullspace(Matrix(transformation.matrix)))' * InternalBoundingBoxes.half_size(bounding_box)))

function random_surface_positions(transformation::Shift{N}, density::GeometryProperties,
    bounding_box::InternalBoundingBox{N}, scale_density) where {N}
    positions, intersections = random_surface_positions(transformation.geometry, density,
        forward(transformation, bounding_box), scale_density)
    [backward(transformation, position) for position in positions], intersections
end

function random_surface_positions(transformation::Scale{N}, density::GeometryProperties,
    bounding_box::InternalBoundingBox{N}, scale_density) where {N}
    positions, intersections = random_surface_positions(transformation.geometry, density,
        forward(transformation, bounding_box), scale_density * transformation.scale^(N - 1))
    [backward(transformation, position) for position in positions], intersections
end

function random_surface_positions(transformation::Rotate{N, M}, density::GeometryProperties,
    bounding_box::InternalBoundingBox{N}, scale_density) where {N, M}
    positions, intersections = random_surface_positions(transformation.geometry, density,
        forward(transformation, bounding_box), scale_density * (N == M ? 1 : _projected_scale(transformation, bounding_box)))
    positions = N == M ? [backward(transformation, position) for position in positions] : _deproject_positions(transformation, positions, bounding_box)
    intersections = [Intersection(hit.distance, backward_normal(transformation, hit.normal), hit.inside,
        hit.obstruction_index, hit.hit_gap) for hit in intersections]
    positions, intersections
end

_geometry_mesh(transformation::Shift, geometry; kwargs...) = _geometry_mesh_preserving(transformation, geometry; kwargs...)
_geometry_mesh(transformation::Scale, geometry; kwargs...) = _geometry_mesh_preserving(transformation, geometry; kwargs...)
_geometry_mesh(transformation::Rotate{N, N}, geometry; kwargs...) where N = _geometry_mesh_preserving(transformation, geometry; kwargs...)
_geometry_mesh(transformation::Shift; kwargs...) = _geometry_mesh(transformation, transformation.geometry; kwargs...)
_geometry_mesh(transformation::Scale; kwargs...) = _geometry_mesh(transformation, transformation.geometry; kwargs...)
_geometry_mesh(transformation::Rotate; kwargs...) = _geometry_mesh(transformation, transformation.geometry; kwargs...)

function _geometry_mesh_preserving(transformation, geometry; bounding_box=nothing, kwargs...)
    local_box = isnothing(bounding_box) ? nothing : forward(transformation, bounding_box)
    _transform_native(transformation, _geometry_mesh(geometry; kwargs..., bounding_box=local_box))
end

_transform_native(transformation, value) = value isa Real ? backward(transformation, SVector{1, Float64}(value))[1] :
    value isa SVector ? backward(transformation, value) :
    value isa NamedTuple ? _mesh_result([backward(transformation, vertex) for vertex in value.vertices], value.triangles) :
    [_transform_native(transformation, item) for item in value]

function _geometry_mesh(transformation::Rotate{3, 1}, geometry; height=nothing, bounding_box=nothing, kwargs...)
    local_box = isnothing(bounding_box) ? nothing : forward(transformation, bounding_box)
    positions = _geometry_mesh(geometry; kwargs..., bounding_box=local_box)
    normal = SVector{3, Float64}(transformation.matrix[1, :])
    reference = abs(normal[1]) < 0.9 ? SVector(1., 0., 0.) : SVector(0., 1., 0.)
    first = cross(normal, reference); first = first ./ norm(first)
    second = cross(normal, first)
    mesh_height(axis) = isnothing(height) ? (isnothing(bounding_box) ? 1. :
        2 * sum(abs.(axis) .* InternalBoundingBoxes.half_size(bounding_box))) : Float64(height)
    h1 = mesh_height(first)
    h2 = mesh_height(second)
    [_mesh_result([transformation.matrix' * SVector{1, Float64}(position) + a * h1/2 * first + b * h2/2 * second
        for (a, b) in ((1,1),(-1,1),(1,-1),(-1,-1))], [SVector(1,2,3), SVector(4,3,2)]) for position in positions]
end

function _circle_triangles(circle)
    triangles = SVector{3, Int}[]
    for index in 1:length(circle)
        next = mod1(index + 1, length(circle))
        push!(triangles, SVector(index, next, length(circle) + index))
        push!(triangles, SVector(next, length(circle) + next, length(circle) + index))
    end
    triangles
end

function _geometry_mesh(transformation::Rotate{3, 2}, geometry; height=nothing, bounding_box=nothing, kwargs...)
    local_box = isnothing(bounding_box) ? nothing : forward(transformation, bounding_box)
    circles = _geometry_mesh(geometry; kwargs..., bounding_box=local_box)
    first = SVector{3, Float64}(transformation.matrix[1, :]); second = SVector{3, Float64}(transformation.matrix[2, :])
    axis = cross(first, second); axis = axis ./ norm(axis)
    extrusion = isnothing(height) ? (isnothing(bounding_box) ? 1. :
        2 * sum(abs.(axis) .* InternalBoundingBoxes.half_size(bounding_box))) : Float64(height)
    [_mesh_result(vcat([transformation.matrix' * point for point in circle] .+ Ref(axis * extrusion/2),
        [transformation.matrix' * point for point in circle] .- Ref(axis * extrusion/2)),
        _circle_triangles(circle)) for circle in circles]
end

end
