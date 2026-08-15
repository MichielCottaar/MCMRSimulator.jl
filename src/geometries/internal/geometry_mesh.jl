"""Render-mesh extraction for fixed physical geometries."""
module GeometryMeshes

import StaticArrays: SVector
import ..Internal: FixedGeometry, geometry_mesh
import ...BoundingBoxes: BoundingBox
import ..InternalBoundingBoxes
import ..PhysicalGeometries: PhysicalGeometry, Transparent, transparent_geometry, InternalBoundingBox
import ..PhysicalGeometries.Groups: GeometryVectorLike, GeometryTuple
import ..PhysicalGeometries.Transformations: Shift, Scale, Rotate, backward, forward
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

function _shift_mesh(mesh_data, displacement)
    [
        _mesh_result(
            [vertex + displacement for vertex in mesh.vertices],
            mesh.triangles,
        )
        for mesh in mesh_data
    ]
end

function geometry_mesh(transformation::Union{Shift{3}, Scale{3}, Rotate{3, 3}}, geometry)
    _transform_mesh(transformation, geometry_mesh(geometry))
end

function geometry_mesh(
    transformation::Union{Shift{3}, Scale{3}, Rotate{3, 3}},
    geometry,
    bounding_box::InternalBoundingBox{3},
)
    local_bounding_box = forward(transformation, bounding_box)
    _transform_mesh(transformation, geometry_mesh(geometry, local_bounding_box))
end

function geometry_mesh(geometry::GeometryVectorLike)
    reduce(vcat, (geometry_mesh(child) for child in geometry); init=NamedTuple[])
end

function geometry_mesh(geometry::GeometryVectorLike, bounding_box::InternalBoundingBox{3})
    reduce(vcat, (geometry_mesh(child, bounding_box) for child in geometry); init=NamedTuple[])
end

function geometry_mesh(geometry::GeometryTuple)
    reduce(vcat, (geometry_mesh(child) for child in geometry); init=NamedTuple[])
end

function geometry_mesh(geometry::GeometryTuple, bounding_box::InternalBoundingBox{3})
    reduce(vcat, (geometry_mesh(child, bounding_box) for child in geometry); init=NamedTuple[])
end

geometry_mesh(repeat::Repeat) = geometry_mesh(repeat.geometry)

function geometry_mesh(repeat::Repeat, bounding_box::InternalBoundingBox{3})
    child_box = InternalBoundingBox(repeat.geometry)
    lower_index = floor.(Int, (InternalBoundingBoxes.lower(bounding_box) - InternalBoundingBoxes.upper(child_box)) ./ repeat.repeats)
    upper_index = ceil.(Int, (InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(child_box)) ./ repeat.repeats)
    child_mesh = geometry_mesh(repeat.geometry)

    result = NamedTuple[]
    for index in Iterators.product((lower_index[i]:upper_index[i] for i in 1:3)...)
        displacement = SVector{3, Float64}(index) .* repeat.repeats
        shifted_lower = InternalBoundingBoxes.lower(child_box) + displacement
        shifted_upper = InternalBoundingBoxes.upper(child_box) + displacement
        all(shifted_lower .<= InternalBoundingBoxes.upper(bounding_box)) || continue
        all(shifted_upper .>= InternalBoundingBoxes.lower(bounding_box)) || continue
        append!(result, _shift_mesh(child_mesh, displacement))
    end
    result
end

geometry_mesh(wrapper::Transparent) = geometry_mesh(transparent_geometry(wrapper))
geometry_mesh(wrapper::Transparent, bounding_box::InternalBoundingBox{3}) =
    geometry_mesh(transparent_geometry(wrapper), bounding_box)
geometry_mesh(geometry::FixedGeometry) = geometry_mesh(geometry.geometry)
geometry_mesh(geometry::FixedGeometry, bounding_box::InternalBoundingBox{3}) =
    geometry_mesh(geometry.geometry, bounding_box)

geometry_mesh(geometry, bounding_box::BoundingBox) =
    geometry_mesh(geometry, InternalBoundingBox(bounding_box))

end
