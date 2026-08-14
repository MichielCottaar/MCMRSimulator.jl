"""Mesh connectivity and component operations."""
module SplitMesh

import SparseArrays: sparse, SparseMatrixCSC
import ..Obstructions: Mesh, nvolumes

"""
    connectivity_matrix(triangles)

Create a sparse matrix connecting vertices that share a triangle edge.
"""
function connectivity_matrix(triangles)
    t1 = [triangle[1] for triangle in triangles]
    t2 = [triangle[2] for triangle in triangles]
    t3 = [triangle[3] for triangle in triangles]
    nvertices = max(maximum(t1), maximum(t2), maximum(t3))
    sparse(vcat(t1, t2, t3), vcat(t2, t3, t1), true, nvertices, nvertices)
end

"""
    connected_indices(connectivity_matrix)
    connected_indices(triangles)

Return the connected-component index for each vertex.
"""
connected_indices(triangles::AbstractVector) = connected_indices(connectivity_matrix(triangles))

function connected_indices(matrix::SparseMatrixCSC)
    nvertices = size(matrix, 1)
    unassigned = trues(nvertices)
    result = zeros(Int, nvertices)
    component_index = 1

    while any(unassigned)
        start = findfirst(unassigned)
        component = falses(nvertices)
        pending = falses(nvertices)
        component[start] = true
        pending[start] = true

        while any(pending)
            for vertex in findall(pending)
                pending[vertex] = false
                for connected in matrix.rowval[matrix.colptr[vertex]:matrix.colptr[vertex + 1] - 1]
                    if !component[connected]
                        component[connected] = true
                        pending[connected] = true
                    end
                end
            end
        end

        @assert all(!(unassigned[vertex] && !component[vertex]) for vertex in 1:nvertices)
        unassigned[component] .= false
        result[component] .= component_index
        component_index += 1
    end

    result
end

"""
    components(mesh)

Return the connected-component index for each mesh triangle. Explicit component
assignments are cached on the user mesh; otherwise they are inferred from mesh
connectivity.
"""
function components(mesh::Mesh)
    mesh.n_obstructions == 0 && return Int[]
    if isnothing(mesh.components.value)
        vertex_indices = connected_indices(mesh.triangles.value)
        mesh.components.value = [vertex_indices[triangle[1]] for triangle in mesh.triangles.value]
    end
    mesh.components.value
end

nvolumes(mesh::Mesh) = maximum(components(mesh); init=0)

end
