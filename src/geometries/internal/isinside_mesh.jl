"""
Special functions to determine if a particle is inside a mesh.

Functions:
- [`isinside_grid`](@ref)
- [`isinside`](@ref)
"""
module IsInsideMesh

import LinearAlgebra: norm, ⋅
import NearestNeighbors: KDTree, nn
import StaticArrays: SVector
import Statistics: mean
import ..Obstructions.Triangles: FullTriangle, IndexTriangle, normal, detect_intersection_partial
import ..FixedObstructionGroups: FixedMesh, detect_intersection_non_repeating, isinside, rotate_from_global, repeating, detect_intersection_loop
import ..BoundingBoxes: BoundingBox
import ..Reflections: Reflection, empty_reflection

"""
    isinside_grid(mesh)

Determines for the centre of each voxel in `mesh.grid`, whether it lies inside or outside of the mesh.
This function assumes that the mesh has been normalised (i.e., all normals point outwards).
"""
function isinside_grid(mesh::FixedMesh)
    mean_triangles = map(o->mean(mesh.vertices[o.indices]), mesh.obstructions)
    if repeating(mesh)
        sz = mesh.grid.size
        mean_triangles = [@. mod(t + sz/2, sz) - sz/2 for t in mean_triangles]
    end
    tree = KDTree(mean_triangles)
    inside_arr = zeros(Bool, size(mesh.grid.indices))
    for index in Tuple.(eachindex(IndexCartesian(), inside_arr))
        centre = @. (index - 0.5) * mesh.grid.resolution + mesh.grid.lower

        triangle_index = nn(tree, centre)[1]
        new_index = triangle_index
        ntry = 0
        while ~iszero(new_index)
            triangle_index = Int(new_index)
            (new_index, _) = detect_intersection_non_repeating(mesh, mean_triangles[triangle_index], centre, triangle_index, true)
            ntry += 1
            if ntry > 1000
                error("Grid voxel centre $centre falls exactly on the edge between triangles. Please shift the mesh a tiny amount to fix this.")
            end
        end
        inpr = (centre - mean_triangles[triangle_index]) ⋅ normal(FullTriangle(mesh.obstructions[triangle_index], mesh.vertices))
        if iszero(inpr)
            @warn "Grid voxel centre falls exactly on the mesh element. This will lead to erroneous estimations of what is the inside of a grid. You might want to shift the mesh a tiny amount."
        end
        inside_arr[index...] = inpr < 0
    end
    return inside_arr
end

function isinside(mesh::FixedMesh{R, O}, pos::SVector{3}, stuck_to::Reflection=empty_reflection) where {R, O}
    rotated = rotate_from_global(mesh, pos)

    if R
        repeats = mesh.grid.size
        half_repeats = repeats/2
        voxel = @. Int(div(rotated, repeats, RoundNearest))
        normed = rotated .- voxel .* repeats
    else
        normed = rotated
    end

    grid_index = @. Int(div(normed - mesh.grid.lower, mesh.grid.resolution, RoundDown)) + 1
    centre = @. (grid_index - 0.5) * mesh.grid.resolution + mesh.grid.lower
    if any(grid_index .< 1) || any(grid_index .> size(mesh.grid.indices))
        return Int[]
    end

    nhit = 0
    for (index, shift_index) in mesh.grid.indices[grid_index...]
        obstruction = FullTriangle(mesh.obstructions[index], mesh.vertices)
        if iszero(shift_index)
            centre_use = centre
            normed_use = normed
        else
            centre_use = centre .- mesh.grid.shifts[shift_index]
            normed_use = normed .- mesh.grid.shifts[shift_index]
        end

        (new_intersection, partial) = detect_intersection_partial(obstruction, centre_use, normed_use)

        if (new_intersection.distance >= 0) && (new_intersection.distance < 1)
            if partial
                nhit += 1
            else
                nhit += 2
            end
        end
    end
    if xor(~iszero(nhit % 4), mesh.grid.isinside[grid_index...])
        return Int[0]
    else
        return Int[]
    end
end
 
end