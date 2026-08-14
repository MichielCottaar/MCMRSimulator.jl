module Mesh

import StaticArrays: SVector, MVector
import LinearAlgebra: cross, norm, svd, ⋅
import DelaunayTriangulation: triangulate, get_triangles

import ...Indices: ObstructionIndex
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, has_inside, InternalBoundingBox
import ..GridDispatch: detect_intersection_grid
import ...InternalBoundingBoxes
import ...InternalBoundingBoxes: grid_indices
import ..BaseObstructions: FullTriangle, normal

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
    grid_bounding_box::InternalBoundingBox{3}
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

function _mesh_boundary_edges(indices)
    edge_occurrences = Dict{Tuple{Int, Int}, Vector{Tuple{Int, Int}}}()
    for triangle in indices
        for edge in ((triangle[1], triangle[2]), (triangle[2], triangle[3]), (triangle[3], triangle[1]))
            key = edge[1] < edge[2] ? edge : (edge[2], edge[1])
            push!(get!(edge_occurrences, key, Tuple{Int, Int}[]), edge)
        end
    end
    [occurrences[1] for occurrences in values(edge_occurrences) if length(occurrences) == 1]
end

function _mesh_boundary_loops(indices)
    boundary_edges = _mesh_boundary_edges(indices)
    outgoing = Dict{Int, Vector{Tuple{Int, Int}}}()
    for edge in boundary_edges
        push!(get!(outgoing, edge[1], Tuple{Int, Int}[]), edge)
    end
    unvisited = Set(boundary_edges)
    loops = Vector{Vector{Int}}()
    while !isempty(unvisited)
        starting_edge = first(unvisited)
        loop = [starting_edge[1]]
        current_edge = starting_edge
        while true
            delete!(unvisited, current_edge)
            push!(loop, current_edge[2])
            current_edge[2] == starting_edge[1] && break
            candidates = get(outgoing, current_edge[2], Tuple{Int, Int}[])
            next_edge = findfirst(edge -> edge in unvisited, candidates)
            isnothing(next_edge) &&
                throw(ArgumentError("mesh boundary edges do not form closed loops"))
            current_edge = candidates[next_edge]
            current_edge[2] in loop[2:end-1] &&
                throw(ArgumentError("mesh boundary edges form a self-intersecting loop"))
        end
        pop!(loop)
        length(loop) >= 3 || throw(ArgumentError("mesh boundary loop must contain at least three vertices"))
        push!(loops, loop)
    end
    loops
end

function _mesh_gap_indices(vertices, indices)
    gap_indices = SVector{3, Int}[]
    for loop in _mesh_boundary_loops(indices)
        positions = vertices[loop]
        centroid = sum(positions) / length(positions)
        centered_positions = hcat((position - centroid for position in positions)...)
        decomposition = svd(centered_positions)
        plane_basis = decomposition.U[:, 1:2]
        projected_positions = [Tuple(plane_basis' * (position - centroid)) for position in positions]
        boundary = [collect(1:length(loop)); 1]
        triangulation = triangulate(
            projected_positions;
            boundary_nodes=[boundary],
            delete_ghosts=true,
            randomise=false,
        )
        for triangle in get_triangles(triangulation)
            all(triangle .> 0) || continue
            push!(gap_indices, SVector{3, Int}(loop[triangle[1]], loop[triangle[2]], loop[triangle[3]]))
        end
    end
    gap_indices
end

function MeshPart(
    vertices,
    indices;
    grid_resolution=nothing,
)
    fixed_vertices = [SVector{3, Float64}(vertex) for vertex in vertices]
    isempty(indices) && throw(ArgumentError("cannot construct a MeshPart without triangles"))
    fixed_indices = [MVector{3, Int}(triangle) for triangle in _mesh_indices(indices, length(fixed_vertices))]
    make_normals_consistent!(fixed_indices)
    gap_indices = _mesh_gap_indices(fixed_vertices, fixed_indices)
    append!(fixed_indices, MVector{3, Int}.(gap_indices))
    first_index_of_gap = length(fixed_indices) - length(gap_indices) + 1
    if curvature(
        fixed_indices,
        fixed_vertices;
        first_index_of_gap,
        include_gap_triangles=true,
    ) < 0
        for triangle in fixed_indices
            triangle[1], triangle[2] = triangle[2], triangle[1]
        end
    end
    fixed_indices = [SVector{3, Int}(triangle) for triangle in fixed_indices]

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
        grid_box,
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

function detect_intersection(
    mesh::MeshPart,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    previous_hit::Intersection{3}=Intersection{3}(),
)
    previous_index = 0
    if !Base.isempty(previous_hit)
        previous_indices = previous_hit.obstruction_index.indices
        length(previous_indices) == 1 ||
            throw(ArgumentError("a non-empty mesh hit must have one triangle index"))
        previous_index = previous_indices[1]
        1 <= previous_index <= length(mesh.indices) ||
            throw(ArgumentError("previous-hit triangle index is not part of the mesh"))
    end

    detect_intersection_grid(
        mesh.grid_bounding_box,
        mesh.inv_resolution,
        mesh.indices_grid,
        start,
        destination,
        previous_hit,
    ) do triangle_index, previous
        triangle_previous_hit = triangle_index == previous_index ? previous : Intersection{3}()
        intersection = detect_intersection(
            triangle(mesh, triangle_index),
            start,
            destination,
            triangle_previous_hit,
        )
        Base.isempty(intersection) && return intersection
        Intersection(
            intersection.distance,
            intersection.normal,
            intersection.inside,
            ObstructionIndex(SVector(triangle_index)),
            triangle_index >= mesh.first_index_of_gap,
        )
    end
end

function _mesh_neighbours(indices)
    norm_edge(edge) = edge[1] > edge[2] ? (edge[2], edge[1]) : edge
    edges(triangle) = norm_edge.(((triangle[1], triangle[2]), (triangle[2], triangle[3]), (triangle[3], triangle[1])))
    edges_to_triangle = Dict{Tuple{Int, Int}, Int}()
    neighbours = Tuple{Int, Int}[]
    for (index, triangle) in enumerate(indices)
        for edge in edges(triangle)
            if haskey(edges_to_triangle, edge)
                push!(neighbours, (edges_to_triangle[edge], index))
            else
                edges_to_triangle[edge] = index
            end
        end
    end
    neighbours
end

function _mesh_curvature_pair(indices, vertices, index1, index2)
    triangle1 = indices[index1]
    triangle2 = indices[index2]
    positions1 = map(index -> vertices[index], triangle1)
    positions2 = map(index -> vertices[index], triangle2)
    position_offset = @. (
        positions1[1] + positions1[2] + positions1[3] -
        positions2[1] - positions2[2] - positions2[3]
    ) / 3
    normal_offset = normal(positions1...) - normal(positions2...)
    (position_offset ⋅ normal_offset) / norm(position_offset)^2
end

function _mesh_triangle_area(indices, vertices, index)
    triangle = indices[index]
    a, b, c = (vertices[vertex] for vertex in triangle)
    norm(cross(b - a, c - a)) / 2
end

"""Compute the area-weighted mean curvature of a triangle mesh."""
function curvature(
    indices,
    vertices;
    first_index_of_gap=length(indices) + 1,
    include_gap_triangles=false,
)
    selected_indices = include_gap_triangles ? indices : indices[1:(first_index_of_gap - 1)]
    neighbours = _mesh_neighbours(selected_indices)
    isempty(neighbours) && return 0.
    weighted_curvatures = [
        (
            _mesh_triangle_area(selected_indices, vertices, index1) +
            _mesh_triangle_area(selected_indices, vertices, index2),
            _mesh_curvature_pair(selected_indices, vertices, index1, index2),
        )
        for (index1, index2) in neighbours
    ]
    valid = filter(pair -> !isnan(pair[2]), weighted_curvatures)
    isempty(valid) && return 0.
    sum(weight * value for (weight, value) in valid) / sum(weight for (weight, _) in valid)
end

curvature(mesh::MeshPart; include_gap_triangles=false) = curvature(
    mesh.indices,
    mesh.vertices;
    first_index_of_gap=mesh.first_index_of_gap,
    include_gap_triangles,
)


"""
    make_normals_consistent!(triangles)

Adjust the triangles to all point outwards or all point inwards.
Assumes that all the triangles are connected (can be enforced using [`connected_components`](@ref)).
"""
function make_normals_consistent!(triangles::AbstractVector)
    edges(t) = ((t[1], t[2]), (t[2], t[3]), (t[3], t[1]))
    counter = Dict{Tuple{Int, Int}, Int}()
    for triangle in triangles
        for (index1, index2) in edges(triangle)
            edge = index2 > index1 ? (index1, index2) : (index2, index1)
            if edge in keys(counter)
                counter[edge] += 1
            else
                counter[edge] = 1
            end
        end
    end
    triangles_seen = fill(false, length(triangles))
    edges_seen = Set{Tuple{Int, Int}}()
    while ~all(triangles_seen)
        starting_triangle = findfirst(.~triangles_seen)
        triangles_seen[starting_triangle] = true
        nseen = 0
        while sum(triangles_seen) > nseen
            nseen = sum(triangles_seen)
            union!(edges_seen, edges(triangles[starting_triangle]))
            for (i, triangle) in enumerate(triangles)
                if triangles_seen[i]
                    continue
                end
                flip_required = -1
                for (index1, index2) in edges(triangle)
                    norm_edge = index2 > index1 ? (index1, index2) : (index2, index1)
                    if counter[norm_edge] != 2
                        continue
                    end
                    if (index1, index2) in edges_seen
                        @assert flip_required in (-1, 1)
                        flip_required = 1
                    elseif (index2, index1) in edges_seen
                        @assert flip_required in (-1, 0)
                        flip_required = 0
                    end
                end
                if flip_required == -1
                    continue
                elseif flip_required == 1
                    v = triangle[1]
                    triangle[1] = triangle[2]
                    triangle[2] = v
                end
                triangles_seen[i] = true
                union!(edges_seen, edges(triangle))
            end
        end
    end
end

end
