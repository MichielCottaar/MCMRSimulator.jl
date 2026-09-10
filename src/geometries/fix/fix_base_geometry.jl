module FixBaseGeometry

import StaticArrays: SVector
import LinearAlgebra: norm
import ...User.Obstructions: Walls, Cylinders, Spheres, FiniteCylinders, Annuli, Mesh
import ...User.SplitMesh: components
import ...User.Obstructions: isglobal
import ...Internal.PhysicalGeometries.BaseObstructions: InfiniteWall, InfiniteCylinder, Sphere, FiniteCylinder
import ...Internal.PhysicalGeometries: Meshes
import ...Internal.PhysicalGeometries.Transformations: Shift

function _values(field_value, number)
    isglobal(field_value) ? fill(field_value.value, number) : collect(field_value.value)
end

function fix_base_geometry(group::Walls)
    fill(InfiniteWall(), length(group))
end

function fix_base_geometry(group::Cylinders)
    [InfiniteCylinder(radius) for radius in _values(group.radius, length(group))]
end

function fix_base_geometry(group::Spheres)
    radii = _values(group.radius, length(group))
    [Sphere(radius) for radius in radii]
end

function fix_base_geometry(group::FiniteCylinders)
    number = length(group)
    positions = [SVector{3, Float64}(position) for position in _values(group.position, number)]
    radii = _values(group.radius, number)
    connected_to = _values(group.connected_to, number)
    use_spherical_endpoint = _values(group.use_spherical_endpoint, number)
    all((connected_to .>= 0) .& (connected_to .<= number)) ||
        throw(ArgumentError("connected_to indices must be between 0 and the number of finite cylinders"))
    all((connected_to .== 0) .| (connected_to .!= collect(1:number))) ||
        throw(ArgumentError("finite cylinders cannot be connected to themselves"))
    spheres = [
        Shift(Sphere(radii[index]), positions[index])
        for index in 1:number if use_spherical_endpoint[index]
    ]
    cylinders = [
        begin 
            displacement = positions[connected_to[index]] - positions[index]
            norm_displacement = norm(displacement)
            # Separate coincident caps at opposing SWC branches from rounding ties.
            additional_displacement = displacement / norm_displacement * sqrt(eps(norm_displacement))
            FiniteCylinder(
                positions[index] - additional_displacement, positions[connected_to[index]] + additional_displacement,
                radii[index], radii[connected_to[index]];
                caps_are_gaps=false,
            )
        end
        for index in 1:number if connected_to[index] != 0
    ]
    spheres, cylinders
end

function fix_base_geometry(group::Annuli)
    inner = [InfiniteCylinder(radius) for radius in _values(group.inner, length(group))]
    outer = [InfiniteCylinder(radius) for radius in _values(group.outer, length(group))]
    (inner, outer)
end

function annuli_size_scale(group::Annuli)
    inner = _values(group.inner, length(group))
    outer = _values(group.outer, length(group))
    minimum((minimum(inner), minimum(outer), minimum(outer .- inner)))
end

function fix_base_geometry(group::Mesh)
    triangle_components = components(group)
    vertices = group.vertices.value
    triangles = group.triangles.value
    unique_components = unique(triangle_components)
    [
        Meshes.Mesh(
            vertices,
            triangles[triangle_components .== component];
            grid_resolution=group.grid_resolution.value,
        )
        for component in unique_components
    ]
end

end
