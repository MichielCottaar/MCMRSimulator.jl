module FixBaseGeometry

import StaticArrays: SVector

import ...User.Obstructions: Walls, Cylinders, Spheres, Annuli, Mesh
import ...User.SplitMesh: components
import ...User.Obstructions: isglobal
import ...InternalNew.PhysicalGeometries.BaseObstructions: InfiniteWall, InfiniteCylinder, Sphere, OverlappingSphere
import ...InternalNew.PhysicalGeometries: Meshes

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
    group.overlapping.value || return [Sphere(radius) for radius in radii]

    positions = [SVector{3, Float64}(position) for position in _values(group.position, length(group))]
    spheres = [OverlappingSphere(radius) for radius in radii]

    for i in eachindex(spheres), j in 1:(i - 1)
        offset = positions[j] - positions[i]
        if sum(offset .* offset) < (radii[i] + radii[j])^2
            push!(spheres[i].overlaps_with, (offset, radii[j]))
            push!(spheres[j].overlaps_with, (-offset, radii[i]))
        end
    end
    spheres
end

function fix_base_geometry(group::Annuli)
    inner = [InfiniteCylinder(radius) for radius in _values(group.inner, length(group))]
    outer = [InfiniteCylinder(radius) for radius in _values(group.outer, length(group))]
    (inner, outer)
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
