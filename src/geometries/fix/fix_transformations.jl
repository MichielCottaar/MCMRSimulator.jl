module FixTransformations

import LinearAlgebra: I, isapprox
import StaticArrays: SMatrix, SVector

import ...User.Obstructions: Walls, Cylinders, Spheres, Annuli, Mesh
import ...User.Obstructions: isglobal
import ...Internal.PhysicalGeometries: PhysicalGeometry
import ...Internal.PhysicalGeometries.Groups: GeometryVector, GeometryTuple
import ...Internal.PhysicalGeometries.Transformations: Shift, Rotate
import ...Internal.PhysicalGeometries.Repeats: Repeat

function _values(field_value, number)
    isglobal(field_value) ? fill(field_value.value, number) : collect(field_value.value)
end

function _shift_obstructions(group, geometry::Vector{P}) where {N, P<:PhysicalGeometry{N}}
    hasproperty(group, :position) || return geometry
    isglobal(group.position) && return geometry
    positions = _values(group.position, length(geometry))
    shifts = [SVector{N, Float64}(position) for position in positions]
    Shift.(geometry, shifts)
end

function _vector_geometry(group, geometry::Vector{P}) where {N, P<:PhysicalGeometry{N}}
    geometry = _shift_obstructions(group, geometry)
    resolution = group.grid_resolution.value
    if !isnothing(resolution)
        return GeometryVector(geometry; grid=true, grid_resolution=resolution)
    elseif group.use_boundingbox.value
        return GeometryVector(geometry; bounding_box=true)
    end
    GeometryVector(geometry)
end

function _apply_local_transformations(group, geometry::Vector)
    geometry = _vector_geometry(group, geometry)
    if !isnothing(group.repeats.value)
        geometry = Repeat(geometry, group.repeats.value)
    end
    geometry
end

function _rotation_matrix(group)
    group.rotation.value
end

function _apply_rotation(group, geometry)
    dimension = group.type.ndim
    matrix = _rotation_matrix(group)
    identity = SMatrix{3, 3, Float64}(I)
    if dimension == 3 && isapprox(matrix, identity; atol=1e-10, rtol=1e-10)
        return geometry
    end
    Rotate(geometry, matrix)
end

function _apply_global_shift(group, geometry)
    hasproperty(group, :position) || return geometry
    isglobal(group.position) || return geometry
    iszero(group.position.value) && return geometry
    position = SVector{group.type.ndim, Float64}(group.position.value)
    Shift(geometry, position)
end

function fix_transformations(group::Union{Walls, Cylinders, Spheres, Mesh}, geometry::Vector)
    geometry = _apply_local_transformations(group, geometry)
    geometry = _apply_global_shift(group, geometry)
    _apply_rotation(group, geometry)
end

function fix_transformations(group::Annuli, geometries::Tuple)
    transformed = map(geometry -> _apply_local_transformations(group, geometry), geometries)
    geometry = GeometryTuple(transformed)
    geometry = _apply_global_shift(group, geometry)
    _apply_rotation(group, geometry)
end

end
