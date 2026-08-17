module Utils

import LinearAlgebra: norm
import StaticArrays: SVector

"""
    icosahedron(subdivision)

Create an icosahedron.
Based on algorithm from https://danielsieger.com/blog/2021/01/03/generating-platonic-solids.html

Each of the 20 original triangles is sub-divided into subdivision^2 triangles
"""
function icosahedron()
    golden_ratio = (1 + sqrt(5)) / 2
    p = 1 / golden_ratio
    vertices = [v / norm(v) for v in [
        [0, p, -1],
        [p, 1, 0],
        [-p, 1, 0],
        [0, p, 1],
        [0, -p, 1],
        [-1, 0, p],
        [0, -p, -1],
        [1, 0, -p],
        [1, 0, p],
        [-1, 0, -p],
        [p, -1, 0],
        [-p, -1, 0]
    ]]
    triangles = [
        [3, 2, 1],
        [2, 3, 4],
        [6, 5, 4],
        [5, 9, 4],
        [8, 7, 1],
        [7, 10, 1],
        [12, 11, 5],
        [11, 12, 7],
        [10, 6, 3],
        [6, 10, 12],
        [9, 8, 2],
        [8, 9, 11],
        [3, 6, 4],
        [9, 2, 4],
        [10, 3, 1],
        [2, 8, 1],
        [12, 10, 7],
        [8, 11, 7],
        [6, 12, 5],
        [11, 9, 5],
    ]
    return SVector{3, Float64}.(vertices), SVector{3, Int}.(triangles)
end

function icosahedron(subdivisions::Int)
    vertices, base_triangles = icosahedron()
    if isone(subdivisions)
        return vertices, base_triangles
    end

    # place vertices at edges
    edge_map = Dict{Tuple{Int, Int}, Vector{Int}}()
    for triangle in base_triangles
        for (i1, i2) in [
            (triangle[1], triangle[2]),
            (triangle[2], triangle[3]),
            (triangle[3], triangle[1]),
        ]
            edge = (i1, i2)
            if edge in keys(edge_map)
                continue
            end
            edge_map[edge] = length(vertices) .+ (1:(subdivisions - 1))
            edge_map[(i2, i1)] = reverse(edge_map[edge])
            intermediates = [vertices[i1] .* (1 - d) .+ vertices[i2] .* d for d in (0:1/(subdivisions+1):1)[2:end-1]]
            append!(vertices, intermediates)
        end
    end

    triangles = SVector{3, Int}[]
    for triangle in base_triangles
        vertex_indices = Dict{Tuple{Int, Int}, Int}()
        for layer in 1:subdivisions+1
            for index in 1:(subdivisions + 2 - layer)
                if layer == 1
                    if index == 1
                        vertex_indices[(layer, index)] = triangle[1]
                    elseif index == subdivisions + 1
                        vertex_indices[(layer, index)] = triangle[2]
                    else
                        vertex_indices[(layer, index)] = edge_map[(triangle[1], triangle[2])][index - 1]
                    end
                elseif layer == subdivisions + 1
                    vertex_indices[(layer, index)] = triangle[3]
                elseif index == 1
                    vertex_indices[(layer, index)] = edge_map[(triangle[1], triangle[3])][layer - 1]
                elseif layer == subdivisions + 2 - layer
                    vertex_indices[(layer, index)] = edge_map[(triangle[2], triangle[3])][layer - 1]
                else
                    v1 = vertices[edge_map[(triangle[1], triangle[3])][layer - 1]]
                    v2 = vertices[edge_map[(triangle[2], triangle[3])][layer - 1]]
                    dist = (index - 1) / (subdivisions + 1 - layer)
                    push!(vertices, (1 - dist) .* v1 .+ dist .* v2)
                    vertex_indices[(layer, index)] = length(vertices)
                end
            end
        end

        for layer in 1:subdivisions
            for index in 1:(subdivisions + 1 - layer)
                push!(triangles, SVector{3, Int}((
                    vertex_indices[(layer, index)],
                    vertex_indices[(layer, index+1)],
                    vertex_indices[(layer+1, index)],
                )))
            end
        end

        for layer in 1:subdivisions-1
            for index in 2:(subdivisions + 1 - layer)
                push!(triangles, SVector{3, Int}((
                    vertex_indices[(layer, index)],
                    vertex_indices[(layer+1, index)],
                    vertex_indices[(layer+1, index-1)],
                )))
            end
        end
    end

    return [v ./ norm(v) for v in vertices], triangles
end

end
