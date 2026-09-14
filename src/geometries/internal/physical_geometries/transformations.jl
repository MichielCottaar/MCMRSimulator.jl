"""Transformations between physical geometry coordinate systems."""
module Transformations

import StaticArrays: SMatrix, SVector
import LinearAlgebra: norm, nullspace, cross
import Random: rand
import ...InternalBoundingBoxes
import ....BoundingBoxes: BoundingBoxNotSupported
import ..PhysicalGeometries: PhysicalGeometry, child_type, find_intersection, find_intersection_requires_inside, get_child, get_intersection_params_requires_inside, has_inside, has_single_inside, inside_indices_eltype, intersection_type, bound_intersection_type, isinside_single, inside_indices, InternalBoundingBox
import ..PhysicalGeometries: random_surface_positions, size_scale, distance_to_surface, _geometry_mesh, _mesh_result, _translate_native
import ..PhysicalGeometries: IntersectionParams, to_child_coordinates, from_child_coordinates, to_child_coordinates_normal
import ...Properties: GeometryProperties

"""
    Transformation{N, M, P} <: PhysicalGeometry{N}

Transformation that places an `M`-dimensional child geometry in an
`N`-dimensional parent coordinate system.

`to_child_coordinates` maps from the parent coordinate system to the child
coordinate system. `from_child_coordinates` maps child coordinates to the
parent coordinate system.
"""
abstract type Transformation{N, M, P<:PhysicalGeometry{M}} <: PhysicalGeometry{N} end

get_intersection_params_requires_inside(::Type{<:Transformation{N, M, P}}) where {N, M, P} =
    get_intersection_params_requires_inside(P)

function Base.show(io::IO, ::Type{T}) where {T <: Transformation}
    print(io, nameof(T))
    T isa UnionAll && return
    print(io, "{")
    show(io, Base.unwrap_unionall(supertype(T)).parameters[end])
    print(io, "}")
end

child_type(::Type{<:Transformation{N, M, P}}) where {N, M, P} = P
inside_indices_eltype(::Type{<:Transformation{N, M, P}}) where {N, M, P} = inside_indices_eltype(P)
intersection_type(::Type{<:Transformation{N, M, P}}) where {N, M, P} = intersection_type(P)
bound_intersection_type(transformation::Transformation, density) =
    bound_intersection_type(transformation.geometry, density)

has_inside(::Type{<:Transformation{N, M, P}}) where {N, M, P} = has_inside(P)
has_single_inside(::Type{<:Transformation{N, M, P}}) where {N, M, P} = has_single_inside(P)
InternalBoundingBox(transformation::Transformation; kwargs...) =
    from_child_coordinates(
        transformation,
        InternalBoundingBox(transformation.geometry; kwargs...),
    )

function isinside_single(
    transformation::Transformation{N, M},
    position::SVector{N, Float64},
    previous_intersection=nothing,
    ; no_deproject=false,
    kwargs...,
) where {N, M}
    no_deproject && transformation isa Rotate && N > M && !all(
        abs.(nullspace(Matrix(transformation.matrix'))' * position) .<= 0.5,
    ) && return false
    isinside_single(
        transformation.geometry,
        to_child_coordinates(transformation, position),
        previous_intersection,
        ; no_deproject, kwargs...,
    )
end

function inside_indices(
    transformation::Transformation{N, M},
    position::SVector{N, Float64},
    intersection=nothing,
    ; no_deproject=false,
    kwargs...,
) where {N, M}
    no_deproject && transformation isa Rotate && N > M && !all(
        abs.(nullspace(Matrix(transformation.matrix'))' * position) .<= 0.5,
    ) && return Tuple[]
    inside_indices(
        transformation.geometry,
        to_child_coordinates(transformation, position),
        intersection,
        ; no_deproject, kwargs...,
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
struct Rotate{N, M, P<:PhysicalGeometry{M}, L} <: Transformation{N, M, P}
    geometry::P
    matrix::SMatrix{N, M, Float64, L}

    function Rotate{N, M, P}(geometry::P, matrix::SMatrix{N, M, Float64}) where {N, M, P<:PhysicalGeometry{M}}
        _validate_rotate_matrix(matrix, N, M)
        new{N, M, P, N*M}(geometry, matrix)
    end
end

function _validate_rotate_matrix(matrix, N, M)
    M <= N || throw(ArgumentError("a Rotate transformation cannot increase dimension"))
    all(
        isapprox(
            sum(matrix[k, i] * matrix[k, j] for k in 1:N),
            i == j ? 1.0 : 0.0;
            atol=1e-10,
            rtol=1e-10,
        )
        for i in 1:M, j in 1:M
    ) || throw(ArgumentError("Rotate matrix columns must be orthonormal"))
    nothing
end

function Rotate(geometry::P, matrix::AbstractMatrix{<:Real}) where {M, P<:PhysicalGeometry{M}}
    N = size(matrix, 1)
    size(matrix, 2) == M || throw(ArgumentError("rotation matrix columns must match geometry dimension"))
    matrix = SMatrix{N, M, Float64}(matrix)
    Rotate{N, M, P}(geometry, matrix)
end

function InternalBoundingBox(
    transformation::Rotate{N, M}; no_deproject=false, kwargs...
) where {N, M}
    if N == M
        return from_child_coordinates(
            transformation,
            InternalBoundingBox(
                transformation.geometry;
                no_deproject,
                kwargs...,
            ),
        )
    end
    no_deproject || throw(BoundingBoxNotSupported(
        "dimension-reducing Rotate transformations do not have a finite bounding box",
    ))
    child_box = InternalBoundingBox(
        transformation.geometry;
        no_deproject,
        kwargs...,
    )
    child_center = InternalBoundingBoxes.center(child_box)
    child_half_size = InternalBoundingBoxes.half_size(child_box)
    center = transformation.matrix * child_center
    null_basis = nullspace(Matrix(transformation.matrix'))
    null_half_size = 0.5 .* vec(sum(abs.(null_basis), dims=2))
    half_size = max.(
        abs.(transformation.matrix) * child_half_size + null_half_size,
        eps(Float64),
    )
    InternalBoundingBoxes.InternalBoundingBox{N}(half_size, center)
end

"""Transform a value from parent coordinates to child coordinates."""
function to_child_coordinates end

"""Transform a value from child coordinates to parent coordinates."""
function from_child_coordinates end

to_child_coordinates(transformation::Shift, position) = position - transformation.shift
from_child_coordinates(transformation::Shift, position) = position + transformation.shift
from_child_coordinates(::Shift, params::IntersectionParams) = params
to_child_coordinates_normal(::Shift, normal) = normal

to_child_coordinates(transformation::Scale, position) = position / transformation.scale
from_child_coordinates(transformation::Scale, position) = transformation.scale * position
from_child_coordinates(::Scale, params::IntersectionParams) = params
to_child_coordinates_normal(::Scale, normal) = normal

distance_to_surface(
    transformation::Transformation,
    position::SVector,
) = distance_to_surface(
    transformation.geometry,
    to_child_coordinates(transformation, position),
)

distance_to_surface(
    transformation::Scale,
    position::SVector,
) = transformation.scale * distance_to_surface(
    transformation.geometry,
    to_child_coordinates(transformation, position),
)

to_child_coordinates(transformation::Shift, box::InternalBoundingBoxes.InternalBoundingBox) =
    InternalBoundingBoxes.shift(box, -transformation.shift)

from_child_coordinates(transformation::Shift, box::InternalBoundingBoxes.InternalBoundingBox) =
    InternalBoundingBoxes.shift(box, transformation.shift)

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

to_child_coordinates(transformation::Scale, box::InternalBoundingBoxes.InternalBoundingBox) =
    _scale_box(Scale(transformation.geometry, inv(transformation.scale)), box)

from_child_coordinates(transformation::Scale, box::InternalBoundingBoxes.InternalBoundingBox) =
    _scale_box(transformation, box)

to_child_coordinates(transformation::Rotate, position) = transformation.matrix' * position
from_child_coordinates(transformation::Rotate, position) = transformation.matrix * position
to_child_coordinates_normal(transformation::Rotate{N, M}, normal) where {N, M} =
    N == M ? transformation.matrix' * normal : throw(ArgumentError("normal_to_child_coordinates is not defined for dimension-reducing Rotate transformations"))
function from_child_coordinates(
    transformation::Rotate{N, M},
    params::IntersectionParams{M},
) where {N, M}
    normal = transformation.matrix * params.normal
    IntersectionParams{N}(
        params.inside,
        normal ./ norm(normal),
        params.hit_gap,
    )
end

function _rotate_box(transformation::Rotate{N, M}, box::InternalBoundingBoxes.InternalBoundingBox{N}) where {N, M}
    center = transformation.matrix' * _box_center(box)
    half_size = abs.(transformation.matrix') * _box_half_size(box)
    return InternalBoundingBoxes.InternalBoundingBox{M}(half_size, center)
end

to_child_coordinates(transformation::Rotate, box::InternalBoundingBoxes.InternalBoundingBox) =
    _rotate_box(transformation, box)

function from_child_coordinates(
    transformation::Rotate{N, M},
    box::InternalBoundingBoxes.InternalBoundingBox{K},
) where {N, M, K}
    N == M && K == M || throw(ArgumentError("from_child_coordinates is not defined for dimension-reducing Rotate transformations"))
    InternalBoundingBoxes.InternalBoundingBox{N}(
        abs.(transformation.matrix) * _box_half_size(box),
        transformation.matrix * _box_center(box),
    )
end

function find_intersection(
    transformation::Transformation{N, M, P},
    start::SVector{N, Float64},
    destination::SVector{N, Float64},
    previous_hit=nothing,
    inside=nothing,
) where {N, M, P}
    find_intersection(
        transformation.geometry,
        to_child_coordinates(transformation, start),
        to_child_coordinates(transformation, destination),
        previous_hit,
        inside,
    )
end

find_intersection_requires_inside(::Type{<:Transformation{N, M, P}}) where {N, M, P} =
    find_intersection_requires_inside(P)

get_child(transformation::Transformation, indices) = (transformation.geometry, indices)

size_scale(geometry::Shift) = size_scale(geometry.geometry)
size_scale(geometry::Rotate) = size_scale(geometry.geometry)
size_scale(geometry::Scale) = geometry.scale * size_scale(geometry.geometry)

function _deproject_positions(transformation::Rotate{N, M}, positions,
    bounding_box::InternalBoundingBox{N}) where {N, M}
    basis = nullspace(Matrix(transformation.matrix'))
    center = InternalBoundingBoxes.center(bounding_box)
    half_size = InternalBoundingBoxes.half_size(bounding_box)
    null_center = basis' * center
    null_half_size = abs.(basis)' * half_size
    [transformation.matrix * position + basis * (null_center +
        (2 .* rand(length(null_center)) .- 1) .* null_half_size) for position in positions]
end

_projected_scale(transformation::Rotate{N, M}, bounding_box::InternalBoundingBox{N}) where {N, M} =
    prod(2 .* (abs.(nullspace(Matrix(transformation.matrix')))' * InternalBoundingBoxes.half_size(bounding_box)))

_surface_density_scale(::Transformation) = 1.0
_surface_density_scale(transformation::Scale{N}) where {N} = transformation.scale^(N - 1)

function _random_surface_positions_same_dimension(
    transformation::Transformation{N, N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density;
    include_gap=false,
    kwargs...,
) where {N}
    positions, indices = random_surface_positions(
        transformation.geometry,
        density,
        to_child_coordinates(transformation, bounding_box),
        scale_density * _surface_density_scale(transformation),
        ; include_gap, kwargs...,
    )
    [from_child_coordinates(transformation, position) for position in positions], indices
end

random_surface_positions(
    transformation::Transformation{N, N},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density;
    include_gap=false,
    kwargs...,
) where {N} = _random_surface_positions_same_dimension(
    transformation, density, bounding_box, scale_density; include_gap, kwargs...,
)

function random_surface_positions(
    transformation::Rotate{N, M},
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{N},
    scale_density;
    include_gap=false,
    no_deproject=false,
    kwargs...,
) where {N, M}
    N == M && return _random_surface_positions_same_dimension(
        transformation, density, bounding_box, scale_density;
        include_gap, no_deproject, kwargs...,
    )
    if no_deproject
        positions, indices = random_surface_positions(
            transformation.geometry,
            density,
            InternalBoundingBox(transformation.geometry; no_deproject=true, kwargs...),
            scale_density;
            include_gap, kwargs...,
        )
        return [transformation.matrix * position for position in positions], indices
    end
    positions, indices = random_surface_positions(
        transformation.geometry,
        density,
        to_child_coordinates(transformation, bounding_box),
        scale_density * _projected_scale(transformation, bounding_box),
        ; include_gap, kwargs...,
    )
    _deproject_positions(transformation, positions, bounding_box), indices
end

_geometry_mesh(transformation::Shift, geometry; kwargs...) = _geometry_mesh_preserving(transformation, geometry; kwargs...)
_geometry_mesh(transformation::Scale, geometry; kwargs...) = _geometry_mesh_preserving(transformation, geometry; kwargs...)
_geometry_mesh(transformation::Rotate{N, N}, geometry; kwargs...) where N = _geometry_mesh_preserving(transformation, geometry; kwargs...)
_geometry_mesh(transformation::Shift; kwargs...) = _geometry_mesh(transformation, transformation.geometry; kwargs...)
_geometry_mesh(transformation::Scale; kwargs...) = _geometry_mesh(transformation, transformation.geometry; kwargs...)
_geometry_mesh(transformation::Rotate; kwargs...) = _geometry_mesh(transformation, transformation.geometry; kwargs...)

function _geometry_mesh_preserving(transformation, geometry; bounding_box=nothing, kwargs...)
    local_box = isnothing(bounding_box) ? nothing : to_child_coordinates(transformation, bounding_box)
    _transform_native(transformation, _geometry_mesh(geometry; kwargs..., bounding_box=local_box))
end

_transform_native(transformation, value) = value isa Real ? from_child_coordinates(transformation, SVector{1, Float64}(value))[1] :
    value isa SVector ? from_child_coordinates(transformation, value) :
    value isa NamedTuple ? _mesh_result([from_child_coordinates(transformation, vertex) for vertex in value.vertices], value.triangles) :
    [_transform_native(transformation, item) for item in value]

function _geometry_mesh(transformation::Rotate{3, 1}, geometry; height=nothing, bounding_box=nothing, kwargs...)
    local_box = isnothing(bounding_box) ? nothing : to_child_coordinates(transformation, bounding_box)
    positions = _geometry_mesh(geometry; kwargs..., bounding_box=local_box)
    normal = SVector{3, Float64}(transformation.matrix[:, 1])
    reference = abs(normal[1]) < 0.9 ? SVector(1., 0., 0.) : SVector(0., 1., 0.)
    first = cross(normal, reference); first = first ./ norm(first)
    second = cross(normal, first)
    mesh_height(axis) = isnothing(height) ? (isnothing(bounding_box) ? 1. :
        2 * sum(abs.(axis) .* InternalBoundingBoxes.half_size(bounding_box))) : Float64(height)
    h1 = mesh_height(first)
    h2 = mesh_height(second)
    [_mesh_result([transformation.matrix * SVector{1, Float64}(position) + a * h1/2 * first + b * h2/2 * second
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
    local_box = isnothing(bounding_box) ? nothing : to_child_coordinates(transformation, bounding_box)
    circles = _geometry_mesh(geometry; kwargs..., bounding_box=local_box)
    first = SVector{3, Float64}(transformation.matrix[:, 1]); second = SVector{3, Float64}(transformation.matrix[:, 2])
    axis = cross(first, second); axis = axis ./ norm(axis)
    extrusion = isnothing(height) ? (isnothing(bounding_box) ? 1. :
        2 * sum(abs.(axis) .* InternalBoundingBoxes.half_size(bounding_box))) : Float64(height)
    [_mesh_result(vcat([transformation.matrix * point for point in circle] .+ Ref(axis * extrusion/2),
        [transformation.matrix * point for point in circle] .- Ref(axis * extrusion/2)),
        _circle_triangles(circle)) for circle in circles]
end

end
