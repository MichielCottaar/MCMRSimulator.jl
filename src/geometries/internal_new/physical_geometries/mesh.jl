module Mesh

import StaticArrays: SVector

import ..PhysicalGeometries: PhysicalGeometry, has_inside, InternalBoundingBox
import ...InternalBoundingBoxes
import ...InternalBoundingBoxes: grid_indices
import ..BaseObstructions: FullTriangle

"""A connected mesh component backed by lazily materialized triangles.

`indices` contains physical mesh triangles followed by triangles that close a
mesh gap. `first_index_of_gap` marks the first gap triangle. All indices use
one-based indices into `vertices`.
"""
struct MeshPart <: PhysicalGeometry{3}
    vertices::Vector{SVector{3, Float64}}
    indices::Vector{SVector{3, Int}}
    first_index_of_gap::Int
    bounding_box::InternalBoundingBox{3}
    inv_resolution::SVector{3, Float64}
    indices_grid::Array{Vector{Int}, 3}
end

has_inside(::Type{MeshPart}) = true

function _mesh_indices(indices, nvertices)
    result = [SVector{3, Int}(triangle) for triangle in indices]
    all(all(triangle .>= 1) && all(triangle .<= nvertices) for triangle in result) ||
        throw(ArgumentError("mesh triangle indices must refer to existing vertices"))
    result
end

function _mesh_grid_box(box::InternalBoundingBox{3})
    lower_bound = InternalBoundingBoxes.lower(box)
    upper_bound = InternalBoundingBoxes.upper(box)
    extent = upper_bound - lower_bound
    # A planar mesh can have zero extent in one or more dimensions. Give the
    # dispatch grid a finite cell in those dimensions without changing the
    # mesh's reported bounding box.
    grid_extent = max.(extent, eps(Float64))
    InternalBoundingBox(
        grid_extent / 2,
        (lower_bound + upper_bound) / 2,
    )
end

function _mesh_grid_resolution(extent, number, resolution)
    isnothing(resolution) && return maximum(extent) / max(ceil(Int, number^(1 / 3)), 1)
    isinf(resolution) && return resolution
    resolution > 0 || throw(ArgumentError("grid_resolution must be positive"))
    Float64(resolution)
end

function MeshPart(
    vertices,
    indices;
    first_index_of_gap=length(indices) + 1,
    grid_resolution=nothing,
)
    fixed_vertices = [SVector{3, Float64}(vertex) for vertex in vertices]
    isempty(indices) && throw(ArgumentError("cannot construct a MeshPart without triangles"))
    fixed_indices = _mesh_indices(indices, length(fixed_vertices))
    1 <= first_index_of_gap <= length(fixed_indices) + 1 ||
        throw(ArgumentError("first_index_of_gap must identify an index or the position after the final index"))

    triangle_boxes = [
        InternalBoundingBox(FullTriangle(fixed_vertices[triangle[1]], fixed_vertices[triangle[2]], fixed_vertices[triangle[3]]))
        for triangle in fixed_indices
    ]
    lower_bound = reduce((a, b) -> min.(a, b), (InternalBoundingBoxes.lower(box) for box in triangle_boxes))
    upper_bound = reduce((a, b) -> max.(a, b), (InternalBoundingBoxes.upper(box) for box in triangle_boxes))
    bounding_box = InternalBoundingBox(
        (upper_bound - lower_bound) / 2,
        (upper_bound + lower_bound) / 2,
    )

    grid_box = _mesh_grid_box(bounding_box)
    extent = InternalBoundingBoxes.upper(grid_box) - InternalBoundingBoxes.lower(grid_box)
    resolution = _mesh_grid_resolution(extent, length(triangle_boxes), grid_resolution)
    dimensions = isinf(resolution) ? SVector{3, Int}(1, 1, 1) :
        max.(SVector{3, Int}(ceil.(Int, extent ./ resolution)), 1)
    cell_size = extent ./ dimensions
    inv_resolution = SVector{3, Float64}(1 ./ cell_size)

    indices_grid = grid_indices(
        grid_box,
        dimensions,
        triangle_boxes,
    )
    MeshPart(
        fixed_vertices,
        fixed_indices,
        first_index_of_gap,
        bounding_box,
        inv_resolution,
        indices_grid,
    )
end

function _mesh_triangle(mesh::MeshPart, triangle::SVector{3, Int})
    FullTriangle(
        mesh.vertices[triangle[1]],
        mesh.vertices[triangle[2]],
        mesh.vertices[triangle[3]],
    )
end

triangle(mesh::MeshPart, index::Int) = _mesh_triangle(mesh, mesh.indices[index])

InternalBoundingBox(mesh::MeshPart) = mesh.bounding_box

end
