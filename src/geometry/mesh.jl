"""
    Mesh(vertices, triangles)

An [`Obstruction`](@ref) formed from a triangular mesh.
"""
struct Mesh{N} <: Obstruction
    vertices :: Vector{PosVector}
    triangles :: Vector{SVector{3, Int}}
    normals :: Vector{PosVector}
    dist_planes :: Vector{Float}
    shape :: GridShape
    grid :: Array{Vector{Int}, 3}
    function Mesh(vertices, triangles, grid_size=20)
        vertices = map(PosVector, vertices)
        triangles = map(SVector{3, Int}, triangles)
        normals = map(t -> normal(vertices[t[1]], vertices[t[2]], vertices[t[3]]), triangles)
        dist_planes = map((n, t) -> n ⋅ vertices[t[1]], normals, triangles)
        bounding_box = expand(BoundingBox(min.(vertices...), max.(vertices...)), 1.001)
        shape = GridShape(bounding_box, grid_size)
        grid = mesh_grid_intersection(shape, vertices, triangles)
        new{length(triangles)}(vertices, triangles, normals, dist_planes, shape, grid)
    end
end

function box_mesh(;center=SA[0, 0, 0], size=[1, 1, 1], grid_size=10)
    center = PosVector(center)
    size = PosVector(size)
    vertices = [
        [0., 0, 0],
        [1., 0, 0],
        [0., 1, 0],
        [1., 1, 0],
        [0., 0, 1],
        [1., 0, 1],
        [0., 1, 1],
        [1., 1, 1],
    ]
    triangles = [
        # -z layer
        [1, 2, 3],
        [2, 3, 4],
        # +z layer
        [5, 6, 7],
        [6, 7, 8],
        # -x layer
        [1, 7, 5],
        [1, 3, 7],
        # +x layer
        [2, 8, 6],
        [2, 4, 8],
        # -y layer
        [1, 2, 6],
        [1, 5, 6],
        # +y layer
        [3, 4, 8],
        [3, 7, 8],
    ]
    c2 = center .- (size .* 0.5)
    Mesh(map(v->PosVector((v .* size) .+ c2), vertices), triangles, grid_size)
end


GridShape(mesh::Mesh) = mesh.shape
BoundingBox(mesh::Mesh) = BoundingBox(GridShape(mesh))

"""
    normal(p1, p2, p3)

Computes the normal of a triangle formed by the three points
"""
function normal(p1 :: PosVector, p2 :: PosVector, p3 :: PosVector)
    e1 = p2 - p1
    e2 = p3 - p1
    tonorm = cross(e1, e2)
    return tonorm ./ norm(tonorm)
end


"""
    mesh_grid_intersection(shape, vertices, triangles)

Determines for each voxel in the grid with which triangles it intersects.

# Algorithm
We follow the algorithm proposed by https://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/tribox_tam.pdf,
which is based on the separating axis theorem.

We project each voxel in the grid and each triangle on a series of 9 axes:
- the normals of the cube ([1, 0, 0], [0, 1, 0], [0, 0, 1])
- the normal of the triangle
- the cross products of any edge of the cube (same as normals) and the edges of the triangle (total of 9 tests)
If the cube and triangle do not overlap along any of these axes, they are non-overlapping.
If they overlap along all axes, the triangle and cube overlap
"""
function mesh_grid_intersection(shape::GridShape, vertices::Vector{PosVector}, triangles :: Vector{SVector{3, Int}})
    sz = shape.size
    grid = [Int[] for _ in 1:sz[1], _ in 1:sz[2], _ in 1:sz[3]]
    hit = fill(true, sz...)
    vertices_voxel = map(v->project(v, shape), vertices)
    for (store_index, index_triangle) in enumerate(triangles)
        fill!(hit, true)
        triangle = map(idx -> vertices_voxel[idx], index_triangle)

        # Check overlap along grid normals
        selector = Any[Colon(), Colon(), Colon()]
        for dim in 1:3
            (low, high) = extrema((triangle[1][dim], triangle[2][dim], triangle[3][dim]))
            low_range = 1:(Int(floor(low)) - 1)
            upper_range = Int(ceil(high)):sz[dim]
            for r in (low_range, upper_range)
                selector[dim] = low_range
                hit[selector...] .= false
                hit[selector...] .= false
            end
            selector[dim] = Colon()
        end

        dims = (
            (1, 2, 3),
            (2, 3, 1),
            (3, 1, 2),
        )

        # check overlap along the other lines
        other_lines = [normal(triangle...)]
        res = [0., 0., 0.]
        for (_, edge1, edge2) in dims
            edge = triangle[edge2] - triangle[edge1]
            for (dim, d1, d2) in dims
                if iszero(edge[d1]) && iszero(edge[d2])
                    continue
                end
                res[dim] = 0.
                res[d1] = -edge[d2]
                res[d2] = edge[d1]
                push!(other_lines, PosVector(res))
            end
        end

        for project_onto in other_lines
            lower, upper = extrema(map(t -> t ⋅ project_onto, triangle))
            diff = sum(abs.(project_onto))
            to_max = 0.5 * sum(project_onto) + diff / 2.
            for index in findall(hit)
                cube_upper = (Tuple(index) ⋅ project_onto) + to_max
                cube_lower = cube_upper - diff
                if cube_upper < lower || cube_lower > upper
                    hit[index] = false
                end
                if (index == CartesianIndex(69, 37, 100)) && (store_index == 3) && false
                    println(project_onto)
                    println(lower, ", ", upper)
                    println(cube_lower, ", ", cube_upper)
                    println(hit[index])
                end
            end
        end
        for index in findall(hit)
            push!(grid[index], store_index)
        end
    end
    grid
end

function detect_collision(movement::Movement, mesh::Mesh{N}, previous::Union{Nothing, Collision}=nothing) where {N}
    collision = nothing
    checked = fill(false, N)
    within_bounds = false
    for (voxel, _, _, grid_time, _) in ray_grid_intersections(GridShape(mesh), movement.origin, movement.destination)
        if any(voxel .< 1) || any(voxel .> size(mesh.grid))
            if within_bounds
                return nothing  # was within bounds, but left it, so we will never return
            end
            continue
        end
        within_bounds = true
        for to_check in mesh.grid[voxel...]
            if checked[to_check]
                continue
            end
            checked[to_check] = true

            tri_solution = detect_collision(movement, mesh, to_check)
            if (
                !isnothing(tri_solution) && 
                (isnothing(previous) || (previous.obstruction !== mesh) || (previous.index != tri_solution.index)) &&
                (isnothing(collision) || tri_solution.distance < collision.distance)
            )
                collision = tri_solution
            end
        end
        if !isnothing(collision) && collision.distance <= grid_time
            return collision
        end
    end
    return nothing
end

function detect_collision(movement::Movement, mesh::Mesh, to_check::Int)
    indices = mesh.triangles[to_check]
    triangle = (
        mesh.vertices[indices[1]],
        mesh.vertices[indices[2]],
        mesh.vertices[indices[3]],
    )
    normal = mesh.normals[to_check]
    dist_plane = mesh.dist_planes[to_check]

    dist_orig = normal ⋅ movement.origin
    dist_dest = normal ⋅ movement.destination
    if dist_orig ≈ dist_dest
        return nothing
    end
    time = (dist_plane - dist_orig) / (dist_dest - dist_orig)
    if time < 0 || time > 1
        return nothing
    end

    intersect_point = time .* movement.destination .+ (1 .- time) .* movement.origin

    dist_plane = normal ⋅ triangle[1]  # distance from origin along normal
    for (dim, d1) in (
            (1, 2),
            (2, 3),
            (3, 1),
        )
        edge = triangle[d1] - triangle[dim]
        to_point = intersect_point - triangle[dim]
        along_normal = cross(edge, to_point)
        if (along_normal ⋅ normal) < 0
            return nothing  # intersect point is on the wrong side of this edge and hence not in the triangle
        end
    end
    return Collision(time, dist_dest > dist_orig ? -normal : normal, mesh, to_check)
end