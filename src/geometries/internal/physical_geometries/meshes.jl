module Meshes

import StaticArrays: SVector, MVector
import LinearAlgebra: cross, norm, svd, ⋅
import DelaunayTriangulation: triangulate, get_triangles
import NearestNeighbors: KDTree, nn

import ...Indices: ObstructionIndex, add_index
import ..PhysicalGeometries: PhysicalGeometry, Intersection, detect_intersection, get_intersection_params, has_inside, has_single_inside, isinside_single, InternalBoundingBox
import ..PhysicalGeometries: intersection_index_length, inside_index_length, all_equal_inside_depth
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh
import ...Properties: GeometryLeafProperties
import ..Groups
import ..GridDispatch: IntersectionGrid, GridIterator, detect_intersection_grid
import ...InternalBoundingBoxes
import ..BaseObstructions: FullTriangle, normal

"""A connected mesh component backed by lazily materialized triangles.

`indices` contains physical mesh triangles followed by triangles that close a
mesh gap. `first_index_of_gap` marks the first gap triangle. All indices use
one-based indices into `vertices`.
"""
struct Mesh <: Groups.GeometryVectorLike{3, FullTriangle}
    vertices::Vector{SVector{3, Float64}}
    indices::Vector{SVector{3, Int}}
    first_index_of_gap::Int
    bounding_box::InternalBoundingBox{3}
    grid::IntersectionGrid{3}
    inside_mask::BitArray{3}
end

struct MeshGeometryView
    mesh::Mesh
end

Groups.group_geometries(mesh::Mesh) = MeshGeometryView(mesh)

Base.size(view::MeshGeometryView) = size(view.mesh.indices)
Base.axes(view::MeshGeometryView) = axes(view.mesh.indices)
Base.length(view::MeshGeometryView) = length(view.mesh.indices)
Base.firstindex(view::MeshGeometryView) = firstindex(view.mesh.indices)
Base.lastindex(view::MeshGeometryView) = lastindex(view.mesh.indices)
Base.getindex(view::MeshGeometryView, index::Int) = triangle(view.mesh, index)
Base.iterate(view::MeshGeometryView, state...) = begin
    next = isempty(state) ? firstindex(view) : first(state)
    next > lastindex(view) && return nothing
    (triangle(view.mesh, next), (next + 1,))
end
Base.eltype(::Type{MeshGeometryView}) = FullTriangle

has_inside(::Type{Mesh}) = true
has_single_inside(::Type{Mesh}) = true
intersection_index_length(::Type{Mesh}) = 1
inside_index_length(::Type{Mesh}) = 0
all_equal_inside_depth(::Type{Mesh}) = true

function _mesh_indices(indices, nvertices)
    result = [SVector{3, Int}(triangle) for triangle in indices]
    all(all(triangle .>= 1) && all(triangle .<= nvertices) for triangle in result) ||
        throw(ArgumentError("mesh triangle indices must refer to existing vertices"))
    result
end

"Helper function to compute holes in meshes (`_mesh_gap_indices`)"
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

"Helper function to compute holes in meshes (`_mesh_gap_indices`)"
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

"Filles holes in meshes"
function _mesh_gap_indices(vertices, indices)
    gap_indices = SVector{3, Int}[]
    for loop in _mesh_boundary_loops(indices)
        positions = vertices[loop]
        centroid = sum(positions) / length(positions)
        centered_positions = hcat((position - centroid for position in positions)...)
        decomposition = svd(centered_positions)
        plane_basis = decomposition.U[:, 1:2]
        projected_positions = [Tuple(plane_basis' * (position - centroid)) for position in positions]
        signed_area = sum(
            projected_positions[index][1] * projected_positions[mod1(index + 1, length(projected_positions))][2] -
            projected_positions[mod1(index + 1, length(projected_positions))][1] * projected_positions[index][2]
            for index in eachindex(projected_positions)
        )
        if signed_area < 0
            reverse!(loop)
            reverse!(projected_positions)
        end
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

function Mesh(
    vertices,
    indices;
    grid_resolution=nothing,
)
    fixed_vertices = [SVector{3, Float64}(vertex) for vertex in vertices]
    isempty(indices) && throw(ArgumentError("cannot construct a Mesh without triangles"))
    fixed_indices = [MVector{3, Int}(triangle) for triangle in _mesh_indices(indices, length(fixed_vertices))]
    make_normals_consistent!(fixed_indices)
    gap_indices = _mesh_gap_indices(fixed_vertices, fixed_indices)
    append!(fixed_indices, MVector{3, Int}.(gap_indices))
    first_index_of_gap = length(fixed_indices) - length(gap_indices) + 1
    make_normals_consistent!(fixed_indices)
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
    grid = IntersectionGrid(triangle_boxes; resolution=grid_resolution)
    inside_mask = _mesh_inside_mask(
        fixed_vertices,
        fixed_indices,
        grid,
    )
    Mesh(
        fixed_vertices,
        fixed_indices,
        first_index_of_gap,
        grid.bounding_box,
        grid,
        inside_mask,
    )
end

function _mesh_triangle(mesh::Mesh, triangle::SVector{3, Int})
    FullTriangle(
        mesh.vertices[triangle[1]],
        mesh.vertices[triangle[2]],
        mesh.vertices[triangle[3]],
    )
end

triangle(mesh::Mesh, index::Int) = _mesh_triangle(mesh, mesh.indices[index])

function get_intersection_params(
    mesh::Mesh,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    indices::Tuple,
)
    triangle_index = indices[1]
    result = get_intersection_params(
        triangle(mesh, triangle_index),
        start,
        destination,
        indices[2:end],
    )
    (
        inside=result.inside,
        normal=result.normal,
        hit_gap=triangle_index >= mesh.first_index_of_gap,
    )
end

Groups.intersection_candidates(mesh::Mesh, start, destination) =
    (
        (triangle_index, triangle(mesh, triangle_index), dist_all_checked)
        for (triangle_index, dist_all_checked) in GridIterator(mesh.grid, start, destination)
    )

InternalBoundingBox(mesh::Mesh) = mesh.bounding_box

function _mesh_inside_mask(
    vertices,
    indices,
    grid::IntersectionGrid{3},
)
    centroids = [
        sum((vertices[vertex] for vertex in triangle)) / 3
        for triangle in indices
    ]
    tree = KDTree(centroids)
    mask = falses(size(grid.indices))
    lower_bound = InternalBoundingBoxes.lower(grid.grid_bounding_box)
    for coordinate in CartesianIndices(mask)
        grid_coordinate = SVector{3, Float64}(Tuple(coordinate)) .- 0.5
        centre = lower_bound .+ grid_coordinate ./ grid.inv_resolution
        triangle_index = nn(tree, centre)[1]
        destination = centre + 1.001 .* (centroids[triangle_index] - centre)
        intersection = detect_intersection_grid(
            grid,
            centre,
            destination,
            Intersection{3, 1}(),
            Intersection{3, 1}(),
        ) do candidate_index, _
            hit = detect_intersection(
                _mesh_triangle(vertices, indices[candidate_index]),
                centre,
                destination,
            )
            Base.isempty(hit) && return Intersection{3, 1}()
            Intersection(
                hit.distance,
                hit.normal,
                hit.inside,
                add_index(ObstructionIndex(), candidate_index),
                false,
            )
        end
        Base.isempty(intersection) && continue
        hit_index = intersection.obstruction_index.indices[1]
        hit_triangle = _mesh_triangle(vertices, indices[hit_index])
        mask[coordinate] = (centre - centroids[hit_index]) ⋅ normal(hit_triangle) < 0
    end
    mask
end

function _mesh_triangle(vertices, triangle::SVector{3, Int})
    FullTriangle(vertices[triangle[1]], vertices[triangle[2]], vertices[triangle[3]])
end

function _mesh_grid_coordinate(mesh::Mesh, position)
    coordinate = Int.(floor.((position - InternalBoundingBoxes.lower(mesh.grid.grid_bounding_box)) .* mesh.grid.inv_resolution)) .+ 1
    any(coordinate .< 1) || any(coordinate .> size(mesh.grid.indices)) ? nothing : SVector{3, Int}(coordinate)
end

function isinside_single(
    mesh::Mesh,
    position::SVector{3, Float64},
    previous_intersection::Intersection{3}=Intersection{3}(),
)
    if !Base.isempty(previous_intersection)
        previous_indices = previous_intersection.obstruction_index.indices
        length(previous_indices) == 1 ||
            throw(ArgumentError("a non-empty mesh intersection must have one triangle index"))
        1 <= previous_indices[1] <= length(mesh.indices) ||
            throw(ArgumentError("previous-hit triangle index is not part of the mesh"))
        return previous_intersection.inside
    end

    coordinate = _mesh_grid_coordinate(mesh, position)
    isnothing(coordinate) && return false

    lower_bound = InternalBoundingBoxes.lower(mesh.grid.grid_bounding_box)
    centre = lower_bound .+ (SVector{3, Float64}(coordinate) .- 0.5) ./ mesh.grid.inv_resolution
    inside = mesh.inside_mask[coordinate...]
    intersection_points = SVector{3, Float64}[]
    for triangle_index in mesh.grid.indices[coordinate...]
        triangle_intersection = detect_intersection(
            triangle(mesh, triangle_index),
            centre,
            position,
            )
        if !Base.isempty(triangle_intersection) && 0 <= triangle_intersection.distance < 1
            point = centre + triangle_intersection.distance .* (position - centre)
            any(norm(point - existing) < 1e-8 for existing in intersection_points) ||
                push!(intersection_points, point)
        end
    end
    isodd(length(intersection_points)) && (inside = !inside)
    inside
end

function detect_intersection(
    mesh::Mesh,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    previous_hit::Intersection{3, M}=Intersection{3, 1}(),
) where {M}
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
        mesh.grid,
        start,
        destination,
        previous_hit,
        Intersection{3, 1}(),
    ) do triangle_index, previous
        triangle_previous_hit = triangle_index == previous_index ? previous : Intersection{3, 1}()
        intersection = detect_intersection(
            triangle(mesh, triangle_index),
            start,
            destination,
            triangle_previous_hit,
        )
        Base.isempty(intersection) && return Intersection{3, 1}()
        Intersection(
            intersection.distance,
            intersection.normal,
            intersection.inside,
            add_index(ObstructionIndex(), triangle_index),
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

curvature(mesh::Mesh; include_gap_triangles=false) = curvature(
    mesh.indices,
    mesh.vertices;
    first_index_of_gap=mesh.first_index_of_gap,
    include_gap_triangles,
)

size_scale(mesh::Mesh) = begin
    ntriangles = mesh.first_index_of_gap - 1
    ntriangles <= 1 && return Inf
    1 / (2 * curvature(mesh))
end

function random_surface_positions(mesh::Mesh, density::GeometryLeafProperties, bounding_box::InternalBoundingBox{3}, scale_density)
    draws = (let
        values = random_surface_positions(
            triangle(mesh, index), density, bounding_box, scale_density)
        values[1], [add_index(hit, index) for hit in values[2]]
    end for index in 1:(mesh.first_index_of_gap - 1))
    Groups._combine(draws, Val(3))
end

function _geometry_mesh(mesh::Mesh; kwargs...)
    number_of_triangles = mesh.first_index_of_gap - 1
    [(vertices=mesh.vertices, triangles=mesh.indices[1:number_of_triangles])]
end


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
