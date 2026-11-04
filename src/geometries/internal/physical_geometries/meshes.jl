module Meshes

import StaticArrays: SVector, MVector
import LinearAlgebra: cross, norm, svd, ⋅
import Random: rand
import Distributions: Poisson
import DelaunayTriangulation: triangulate, get_triangles

import ..PhysicalGeometries: PhysicalGeometry, IntersectionParams, find_intersection, get_intersection_params, has_inside, has_single_inside, inside_indices_eltype, intersection_type, isinside_single, InternalBoundingBox
import ..PhysicalGeometries: random_surface_positions, size_scale, _geometry_mesh
import ..PhysicalGeometries: inside_sampling_per_obstruction
import ..PhysicalGeometries: leaf_bounding_boxes, finalise!, simplify
import ..PhysicalGeometries: filter_bounding_box, _bounding_boxes_overlap
import ...Properties: GeometryLeafProperties
import ..Groups
import ...InternalBoundingBoxes
import ..BaseObstructions: FullTriangle, normal

"""A connected mesh component backed by lazily materialized triangles.

`indices` contains physical mesh triangles followed by triangles that close a
mesh gap. `first_index_of_gap` marks the first gap triangle. `original_indices`
maps stored triangles to their IDs in the source mesh, or is `nothing` for an
original mesh. `centre_inside` provides a reference point and its inside state.
"""
struct Mesh <: Groups.GeometryVectorLike{3, FullTriangle}
    vertices::Vector{SVector{3, Float64}}
    indices::Vector{SVector{3, Int}}
    first_index_of_gap::Int
    original_indices::Union{Nothing, Vector{Int}}
    centre_inside::Tuple{SVector{3, Float64}, Bool}
end

Base.show(io::IO, ::Type{Mesh}) = print(io, "Mesh")
Base.show(io::IO, ::Mesh) = print(io, "Mesh(...)")

# Mesh triangles are exposed through a generated view, not stored child geometries.
finalise!(::Mesh) = nothing
simplify(mesh::Mesh, ::Val) = mesh

struct MeshGeometryView
    mesh::Mesh
    include_gap::Bool
end

Groups.group_geometries(mesh::Mesh; include_gap=true) = MeshGeometryView(mesh, include_gap)
Groups.indices(mesh::Mesh) = isnothing(mesh.original_indices) ? Base.OneTo(length(mesh.indices)) : mesh.original_indices

leaf_bounding_boxes(mesh::Mesh; inside=false) = inside ?
    InternalBoundingBox(mesh) :
    [leaf_bounding_boxes(_triangle_at(mesh, index)) for index in eachindex(mesh.indices)]

Base.size(view::MeshGeometryView) = (length(view),)
Base.axes(view::MeshGeometryView) = (Base.OneTo(length(view)),)
Base.length(view::MeshGeometryView) = view.include_gap ? length(view.mesh.indices) : view.mesh.first_index_of_gap - 1
Base.firstindex(view::MeshGeometryView) = firstindex(view.mesh.indices)
Base.lastindex(view::MeshGeometryView) = length(view)
Base.getindex(view::MeshGeometryView, index::Int) = _triangle_at(view.mesh, index)
Base.iterate(view::MeshGeometryView, state=firstindex(view)) = begin
    next = state
    next > lastindex(view) && return nothing
    (_triangle_at(view.mesh, next), next + 1)
end
Base.eltype(::Type{MeshGeometryView}) = FullTriangle

has_inside(::Type{Mesh}) = true
has_single_inside(::Type{Mesh}) = true
inside_indices_eltype(::Type{Mesh}) = Tuple{}
intersection_type(::Type{Mesh}) = Tuple{Int}
function _mesh_indices(indices, nvertices)
    result = [SVector{3, Int}(triangle) for triangle in indices]
    all(all(triangle .>= 1) && all(triangle .<= nvertices) for triangle in result) ||
        throw(ArgumentError("mesh triangle indices must refer to existing vertices"))
    result
end

function _mesh_bounding_box(vertices, empty_center=nothing)
    if isempty(vertices)
        isnothing(empty_center) && throw(ArgumentError("cannot construct a Mesh without vertices"))
        return InternalBoundingBox{3}(zero(SVector{3, Float64}), empty_center)
    end
    lower = reduce((a, b) -> min.(a, b), vertices)
    upper = reduce((a, b) -> max.(a, b), vertices)
    InternalBoundingBox{3}((upper - lower) / 2, (upper + lower) / 2)
end

function _triangle_at(mesh::Mesh, position::Int)
    _mesh_triangle(mesh.vertices, mesh.indices[position])
end

triangle(mesh::Mesh, original_index::Int) = _triangle_at(mesh, Groups._stored_index(mesh, original_index))

function _mesh_subset(mesh::Mesh, selected_indices, centre_inside::Tuple{SVector{3, Float64}, Bool})
    # Keep source triangle IDs while compacting both the triangle list and vertices.
    positions = sort!(unique!([Groups._stored_index(mesh, index) for index in selected_indices]))
    isempty(positions) && return Mesh(SVector{3, Float64}[], SVector{3, Int}[], 1, Int[], centre_inside)

    physical_positions = filter(position -> position < mesh.first_index_of_gap, positions)
    gap_positions = filter(position -> position >= mesh.first_index_of_gap, positions)
    ordered_positions = [physical_positions; gap_positions]
    original_indices = [Groups.indices(mesh)[position] for position in ordered_positions]
    selected_triangles = mesh.indices[ordered_positions]
    vertex_ids = unique!(reduce(vcat, (collect(triangle) for triangle in selected_triangles)))
    vertex_position = Dict(vertex_id => position for (position, vertex_id) in enumerate(vertex_ids))
    vertices = mesh.vertices[vertex_ids]
    triangles = [
        SVector{3, Int}(vertex_position[triangle[1]], vertex_position[triangle[2]], vertex_position[triangle[3]])
        for triangle in selected_triangles
    ]
    Mesh(vertices, triangles, length(physical_positions) + 1, original_indices, centre_inside)
end

function _inside_between(vertices, triangles, start, destination, initial_inside)
    intersections = SVector{3, Float64}[]
    for triangle_indices in triangles
        hit = find_intersection(_mesh_triangle(vertices, triangle_indices), start, destination)
        isnothing(hit) && continue
        _, distance = hit
        0 <= distance < 1 || continue
        point = start + distance .* (destination - start)
        any(norm(point - existing) < 1e-8 for existing in intersections) ||
            push!(intersections, point)
    end
    isodd(length(intersections)) ? !initial_inside : initial_inside
end

function _centre_inside(vertices, triangles)
    bounding_box = _mesh_bounding_box(vertices)
    centre = InternalBoundingBoxes.center(bounding_box)
    upper = InternalBoundingBoxes.upper(bounding_box)
    extent = InternalBoundingBoxes.upper(bounding_box) - InternalBoundingBoxes.lower(bounding_box)
    outside = upper .+ max(maximum(extent), 1.0) .* SVector(1.317, 1.713, 2.117)
    (centre, _inside_between(vertices, triangles, outside, centre, false))
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

"""
    Mesh(vertices, triangle_indices)

Construct a complete mesh component and its gap-closing triangles. Meshes do not
build a private spatial grid: their triangles are traversed directly, and inside
queries start from the bounding-box centre and its computed inside state.
"""
function Mesh(vertices, indices)
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

    Mesh(
        fixed_vertices,
        fixed_indices,
        first_index_of_gap,
        nothing,
        _centre_inside(fixed_vertices, fixed_indices),
    )
end

function get_intersection_params(
    mesh::Mesh,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    indices::Tuple,
    isinside=nothing,
)
    triangle_index = indices[1]
    triangle_position = Groups._stored_index(mesh, triangle_index)
    result = get_intersection_params(
        triangle(mesh, triangle_index),
        start,
        destination,
        indices[2:end],
        isinside,
    )
    IntersectionParams{3}(
        result.inside,
        result.normal,
        triangle_position >= mesh.first_index_of_gap,
    )
end

Groups.intersection_candidates(mesh::Mesh, start, destination) =
    ((index, child, 0.0) for (index, child) in Groups._indexed_children(mesh))

InternalBoundingBox(mesh::Mesh) = _mesh_bounding_box(mesh.vertices, mesh.centre_inside[1])

function _filter_mesh_bounding_box(mesh::Mesh, box::InternalBoundingBox{3})
    centre = InternalBoundingBoxes.center(box)
    centre_inside = isinside_single(mesh, centre)
    selected_indices = Int[]
    for (position, original_index) in enumerate(Groups.indices(mesh))
        triangle_box = InternalBoundingBox(_triangle_at(mesh, position))
        _bounding_boxes_overlap(triangle_box, box) && push!(selected_indices, original_index)
    end
    _mesh_subset(mesh, selected_indices, (centre, centre_inside))
end

function filter_bounding_box(mesh::Mesh, box::InternalBoundingBox{3}, ::Val{false})
    _bounding_boxes_overlap(InternalBoundingBox(mesh), box) || return nothing
    _filter_mesh_bounding_box(mesh, box)
end

function filter_bounding_box(mesh::Mesh, box::InternalBoundingBox{3}, ::Val{true})
    _filter_mesh_bounding_box(mesh, box)
end

function inside_sampling_per_obstruction(
    mesh::Mesh,
    density::Real,
    rng;
    no_deproject=false,
    bounding_box=nothing,
)
    box = InternalBoundingBox(mesh)
    lower_bound = InternalBoundingBoxes.lower(box)
    size = InternalBoundingBoxes.upper(box) - lower_bound
    number_of_samples = rand(rng, Poisson(density * prod(size)))
    positions = SVector{3, Float64}[]
    sizehint!(positions, number_of_samples)
    for _ in 1:number_of_samples
        position = lower_bound + SVector{3, Float64}(rand(rng, 3)) .* size
        isinside_single(mesh, position; no_deproject) && push!(positions, position)
    end
    positions
end

function _mesh_triangle(vertices, triangle::SVector{3, Int})
    FullTriangle(vertices[triangle[1]], vertices[triangle[2]], vertices[triangle[3]])
end

function isinside_single(
    mesh::Mesh,
    position::SVector{3, Float64},
    previous_intersection=nothing,
    ; kwargs...,
)
    if !isnothing(previous_intersection)
        return previous_intersection[2]
    end

    centre, centre_inside = mesh.centre_inside
    _inside_between(mesh.vertices, mesh.indices, centre, position, centre_inside)
end

function find_intersection(
    mesh::Mesh,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    previous_hit=nothing,
    inside=nothing,
    candidate::Union{Nothing, Tuple}=nothing,
)
    previous_index = isnothing(previous_hit) ? nothing : previous_hit[1]
    if !isnothing(candidate)
        length(candidate) == 1 || throw(ArgumentError("a mesh candidate must contain one triangle index"))
        triangle_index = first(candidate)
        child = triangle(mesh, triangle_index)
        triangle_index == previous_index && return nothing
        hit = find_intersection(child, start, destination, nothing, nothing, ())
        return isnothing(hit) ? nothing : (triangle_index, hit...)
    end
    best = nothing
    for (triangle_index, child) in Groups._indexed_children(mesh)
        triangle_index == previous_index && continue
        hit = find_intersection(child, start, destination)
        isnothing(hit) && continue
        inside, distance = hit
        (isnothing(best) || distance < best[3]) || continue
        best = (triangle_index, inside, distance)
    end
    best
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
