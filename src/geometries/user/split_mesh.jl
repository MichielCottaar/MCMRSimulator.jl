"""
Mesh-specific operations to normalise and split a mesh.
- [`split_mesh`](@ref): applies all the normalisations and splitting
"""
module SplitMesh

import StaticArrays: SVector, MVector
import SparseArrays: sparse, SparseMatrixCSC
import LinearAlgebra: norm, ⋅
import Statistics: mean
import ...Internal.Obstructions.Triangles: normal
import ..Obstructions: Mesh, isglobal

"""
    split_mesh(mesh::Mesh)

Returns a vector of meshes that have been adjusted in the following ways:
1. the complete mesh has been split into its connected components.
2. the order of the indices in the triangles has been fixed, so that all normals point outwards.
"""
function split_mesh(old_mesh::Mesh)
    result = Mesh[]
    for (triangle_indices, vertex_indices, triangles) in connected_components(old_mesh.triangles.value)
        make_normals_consistent!(triangles)
        if mean(sign.(curvature(triangles, old_mesh.vertices.value[vertex_indices]))) < 0
            # flip all triangles
            for t in triangles
                v = t[1]
                t[1] = t[2]
                t[2] = v
            end
        end
        kwargs = Dict{Symbol, Any}()
        for key in old_mesh.unique_keys
            if key == :triangles
                kwargs[key] = triangles
            elseif key == :vertices
                kwargs[key] = old_mesh.vertices.value[vertex_indices]
            else
                fv = getproperty(old_mesh, key)
                if isglobal(fv)
                    kwargs[key] = fv.value
                else
                    kwargs[key] = fv.value[triangle_indices]
                end
            end
        end

        push!(result, Mesh(; kwargs...))
    end
    return result
end


"""
    make_normals_consistent!(triangles)

Adjust the triangles to all point outwards or all point inwards.
Assumes that all the triangles are connected (can be enforced using [`connected_components`](@ref)).
"""
function make_normals_consistent!(triangles::AbstractVector)
    edges(t) = [(t[1], t[2]), (t[2], t[3]), (t[3], t[1])]
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

"""
    connectivity_matrix(triangles, nvertices)

Create a nverticesxnvertices sparse boolean matrix, which is true for any connected vertices.
"""
function connectivity_matrix(triangles, nvertices=nothing)
    t1 = [t[1] for t in triangles]
    t2 = [t[2] for t in triangles]
    t3 = [t[3] for t in triangles]
    if isnothing(nvertices)
        nvertices = max(maximum(t1), maximum(t2), maximum(t3))
    end
    sparse(vcat(t1, t2, t3), vcat(t2, t3, t1), true, nvertices, nvertices)
end

"""
    connected_components(triangles)

Splits the mesh represented by `triangles` into a sequences of meshes that are actually internally connected.
"""
function connected_components(triangles::AbstractVector, nvertices=nothing)
    m = connectivity_matrix(triangles, nvertices)
    vertex_indices = connected_components(m)
    indices = [vertex_indices[t[1]] for t in triangles]
    function get_index(i)
        vertex_index = (1:length(vertex_indices))[vertex_indices .== i]
        triangle_index = (1:length(indices))[indices .== i]
        new_index = zeros(Int, length(vertex_indices))
        new_index[vertex_indices .== i] = 1:length(vertex_index)
        return (triangle_index, vertex_index, [MVector{3, Int}(new_index[t]) for t in triangles])
    end
    return get_index.(1:maximum(indices))
end

"""
    connected_components(connectivity_matrix)

Returns a vector of indices, where each connected component has been given a unique index (starting with 1).
"""
function connected_components(m::SparseMatrixCSC)
    nvertices = size(m, 1)
    not_assigned = fill(true, nvertices)
    result = fill(0, nvertices)
    component_index = 1
    while any(not_assigned)
        start_index :: Int = findfirst(not_assigned)
        new_component = BitVector(fill(false, nvertices))
        function mark_index(index)
            if new_component[index]
                return
            end
            new_component[index] = true
            mark_index.(m.rowval[m.colptr[index]:m.colptr[index + 1] - 1])
            nothing
        end
        mark_index(start_index)

        @assert sum(not_assigned .& new_component) == sum(new_component) "new component includes elements already assigned to previous components"
        not_assigned[new_component] .= false
        result[new_component] .= component_index
        component_index += 1
    end
    #@assert all(result .> 0)
    return result
end

"""
    neighbours(triangles)

Return pairs of indices of triangles that share an edge.
"""
function neighbours(triangles::AbstractVector)
    norm_edge((i1, i2)) = i1 > i2 ? (i2, i1) : (i1, i2)
    edges(t) = norm_edge.([(t[1], t[2]), (t[2], t[3]), (t[3], t[1])])

    edges_to_triangle = Dict{Tuple{Int, Int}, Int}()
    pairs = Tuple{Int, Int}[]
    for (i, t) in enumerate(triangles)
        for e in edges(t)
            if e in keys(edges_to_triangle)
                push!(pairs, (edges_to_triangle[e], i))
            else
                edges_to_triangle[e] = i
            end
        end
    end
    return pairs
end

"""
    curvature(triangles, vertices[index1, index2])

Computes the curvature between two neighbouring triangles (`index1` and `index2`).
If no indices are provided computes the curvature between any neighbouring triangles.
"""
function curvature(triangles, vertices, index1, index2)
    pos1 = map(i->vertices[i], triangles[index1])
    pos2 = map(i->vertices[i], triangles[index2])
    pos_offset = @. (pos1[1] + pos1[2] + pos1[3] - pos2[1] - pos2[2] - pos2[3]) / 3
    normal_offset = normal(pos1...) - normal(pos2...)
    return (pos_offset ⋅ normal_offset) / norm(pos_offset)^2
end

curvature(triangles, vertices) = [curvature(triangles, vertices, i1, i2) for (i1, i2) in neighbours(triangles)]


end