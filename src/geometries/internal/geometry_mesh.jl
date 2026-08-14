"""Render-mesh extraction for fixed physical geometries."""
module GeometryMeshes

import ..Internal: FixedGeometry, geometry_mesh
import ..PhysicalGeometries: PhysicalGeometry, Transparent, transparent_geometry
import ..PhysicalGeometries.Groups: GeometryVectorLike, GeometryTuple
import ..PhysicalGeometries.Transformations: Shift, Scale, Rotate, backward
import ..PhysicalGeometries.Repeats: Repeat
import ..PhysicalGeometries.Meshes: Mesh

function _mesh_result(vertices, triangles)
    (vertices=vertices, triangles=triangles)
end

function geometry_mesh(mesh::Mesh)
    number_of_triangles = mesh.first_index_of_gap - 1
    [_mesh_result(mesh.vertices, mesh.indices[1:number_of_triangles])]
end

function _transform_mesh(transformation, mesh_data)
    [
        _mesh_result(
            [backward(transformation, vertex) for vertex in mesh.vertices],
            mesh.triangles,
        )
        for mesh in mesh_data
    ]
end

function geometry_mesh(transformation::Union{Shift{3}, Scale{3}, Rotate{3, 3}}, geometry)
    _transform_mesh(transformation, geometry_mesh(geometry))
end

function geometry_mesh(geometry::GeometryVectorLike)
    reduce(vcat, (geometry_mesh(child) for child in geometry); init=NamedTuple[])
end

function geometry_mesh(geometry::GeometryTuple)
    reduce(vcat, (geometry_mesh(child) for child in geometry); init=NamedTuple[])
end

geometry_mesh(repeat::Repeat) = geometry_mesh(repeat.geometry)
geometry_mesh(wrapper::Transparent) = geometry_mesh(transparent_geometry(wrapper))
geometry_mesh(geometry::FixedGeometry) = geometry_mesh(geometry.geometry)

end
