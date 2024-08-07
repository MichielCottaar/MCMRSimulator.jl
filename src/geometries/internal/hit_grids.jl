module HitGrids

import StaticArrays: SVector
import Statistics: mean
import NearestNeighbors: KDTree, nn
import LinearAlgebra: ⋅
import ..BoundingBoxes: BoundingBox, lower, upper, could_intersect
import ..Obstructions: FixedObstruction, detect_intersection, ObstructionIntersection, isinside, empty_obstruction_intersections
import ..Obstructions.Triangles: normal, detect_intersection_partial, IndexTriangle, FullTriangle
import ..RayGridIntersection: ray_grid_intersections


"""
    HitGrid(obstructions, vertices, grid_resolution[, repeats])

A precomputed grid that makes it easier to compute the interactions between a trajectory and obstructions.

Predetermines the interactions between the obstructions and the grid.
This is determined purely based on the bounding boxes (could be improved by considering the shape of each obstruction).
"""
abstract type HitGrid{N, O} end

BoundingBox(g::HitGrid) = g.bounding_box
inv_resolution(g::HitGrid) = g.inv_resolution
lower(g::HitGrid) = lower(BoundingBox(g))
upper(g::HitGrid) = upper(BoundingBox(g))


"""
    obstructions(grid/fixed_obstructions)

Returns a sequence of all the obstructions included in this group.

This operation is very slow and should be avoided within inner loops.
"""
function obstructions(g::HitGrid{N, O}) where {N, O}
    individual = Dict{Int32, O}()
    for index_list in g.indices
        for pack in index_list
            individual[pack[1]] = pack[end]
        end
    end
    @assert all(sort([keys(individual)...]) .== 1:length(individual))
    return map(key -> individual[key], 1:length(individual))
end


"""
A specific version of [`HitGrid`](@ref) that is optimised for any non-repeating geometry.
"""
struct HitGridNoRepeat{N, O} <: HitGrid{N, O}
    bounding_box :: BoundingBox{N}
    inv_resolution :: SVector{N, Float64}
    indices :: Array{Vector{Tuple{Int32, O}}, N}
end

"""
A specific version of [`HitGrid`](@ref) that is optimised for any repeating geometry.
"""
struct HitGridRepeat{N, O} <: HitGrid{N, O}
    bounding_box :: BoundingBox{N}
    inv_resolution :: SVector{N, Float64}
    indices :: Array{Vector{Tuple{Int32, Int32, O}}, N}
    shifts :: Vector{SVector{N, Float64}}
end

function (::Type{HitGrid})(obstructions::Vector{<:FixedObstruction{N}}, grid_resolution::Float64, repeats::AbstractArray, args...) where {N}
    return HitGrid(obstructions, grid_resolution, SVector{N, Float64}(repeats), args...)
end

function (::Type{HitGrid})(obstructions::Vector{<:FixedObstruction{N}}, grid_resolution::Float64, repeats::Union{Nothing, SVector{N, Float64}}, args...) where {N}
    bounding_boxes = map(o->BoundingBox(o, args...), obstructions)
    bb_actual = BoundingBox(bounding_boxes)
    if ~isnothing(repeats)
        new_lower = map(bb_actual.lower, bb_actual.upper, repeats) do l, u, r
            if u > r/2
                return -r/2
            else
                return l
            end
        end
        new_upper = map(bb_actual.lower, bb_actual.upper, repeats) do l, u, r
            if l < -r/2
                return r/2
            else
                return u
            end
        end
        bb_actual = BoundingBox(new_lower, new_upper)
    end
    extend_by = isfinite(grid_resolution) ? grid_resolution / 100 : 1e-3
    bb = BoundingBox(bb_actual.lower .- extend_by, bb_actual.upper .+ extend_by)
    sz = upper(bb) - lower(bb)
    if isinf(grid_resolution)
        dims = SVector{3, Int}(fill(1, N))
    else
        dims = Int.(div.(sz, grid_resolution, RoundUp))
    end
    actual_grid_resolution = sz ./ dims

    shifts = SVector{N, Float64}[]
    grid = [Tuple{Int32, Int32}[] for _ in Iterators.product(UnitRange.(1, dims)...)]
    for (index, this_bb) in enumerate(bounding_boxes)
        if isnothing(repeats)
            range_repeats = SVector{N, UnitRange{Int64}}(fill(0:0, N))
        else
            range_repeats = UnitRange.(Int.(div.(this_bb.lower, repeats, RoundNearest)), Int.(div.(this_bb.upper, repeats, RoundNearest)))
        end
        for int_shifts in Iterators.product(range_repeats...)
            if isnothing(repeats)
                shift = zero(SVector{N, Float64})
            else
                shift = -SVector{N, Float64}(int_shifts .* repeats)
            end

            l = max.(Int.(div.((lower(this_bb) .+ shift .- lower(bb)), actual_grid_resolution, RoundUp)), 1)
            u = min.(Int.(div.((upper(this_bb) .+ shift .- lower(bb)), actual_grid_resolution, RoundDown)) .+ 1, dims)

            @assert all(u .>= l)

            # found possible intersections
            if all(iszero.(shift))
                index_shift = 0
            else
                if ~(shift in shifts)
                    push!(shifts, shift)
                end
                index_shift = findfirst(s->all(s.==shift), shifts)
            end
            for grid_indices in Iterators.product(UnitRange.(l, u)...)
                push!(grid[grid_indices...], (index, index_shift))
            end
        end
    end
    if length(shifts) == 0
        return HitGridNoRepeat(
            bb,
            dims ./ sz,
            map(grid) do index_arr
                map(i -> (i[1], obstructions[i[1]]), index_arr)
            end,
        )
    else
        return HitGridRepeat(
            bb,
            dims ./ sz,
            map(grid) do index_arr
                map(i -> (i[1], i[2], obstructions[i[1]]), index_arr)
            end,
            shifts,
        )
    end

end


function get_coordinates(grid::HitGrid{N}, position::SVector{N, Float64}) where {N}
    return @. Int(floor((position - lower(grid)) * grid.inv_resolution)) + 1
end

"""
    isinside(grid, position)

Get the indices of obstructions that contain the `position`.
"""
function isinside(grid::HitGrid{N}, position::SVector{N, Float64}, stuck_to::Int32, inside::Bool) where {N}
    res = Int32[]
    for pack in grid.indices[get_coordinates(grid, position)...]
        if grid isa HitGridNoRepeat
            (index, obstruction) = packed
            pos_shift = position
        else
            (index, shift, obstruction) = packed
            if iszero(shift)
                pos_shift = position
            else
                shift = g.grid.shifts[shift_index]
                pos_shift = position .- shift
            end
        end
        if index == stuck_to
            if inside
                push!(res, index)
            end
        elseif isinside(obstruction, pos_shift)
            push!(res, index)
        end
    end
    return res
end


"""
    detect_intersection_grid(grid, start, dest, prev_index, prev_inside, args...)

Computes the first intersection between the trajectory from `start` to `dest` and the obstructions stored in `grid`.

This function returns the index of the obstruction the trajectory intersects with and an [`ObstructionIntersection`](@ref) object describing the intersection in more detail.
If `prev_index` is non-zero, it is assumed that the trajectory starts from hitting the obstruction with that index having hit at `prev_inside`.
"""
function detect_intersection_grid(grid::HitGrid{N}, start::SVector{N, Float64}, dest::SVector{N, Float64}, prev_index::Int32, prev_inside::Bool, args...) where {N}
    if ~could_intersect(BoundingBox(grid), start, dest)
        return (zero(Int32), empty_obstruction_intersections[N])
    end
    if all(size(grid.indices) .== 1)
        return detect_intersection_inner(g, start, dest, prev_index, prev_inside, zero(SVector{N, Int64}) .+ 1, args...)
    end
    found_intersection = empty_obstruction_intersections[N]
    found_intersection_index = zero(Int32)
    valid_voxel = false
    for (voxel, _, _, t2, _) in ray_grid_intersections((start .- lower(grid)) .* grid.inv_resolution, (dest .- lower(grid)) .* grid.inv_resolution)
        if any(voxel .< 0) || any(voxel .>= size(grid.indices))
            if valid_voxel
                break
            else
                continue
            end
        end
        valid_voxel = true
        (intersection_index, intersection) = detect_intersection_inner(grid, start, dest, prev_index, prev_inside, voxel .+ 1, args...)
        if intersection.distance < found_intersection.distance
            found_intersection = intersection
            found_intersection_index = intersection_index
        end
        if found_intersection.distance <= t2
            return (found_intersection_index, found_intersection)
        end
    end
    return (found_intersection_index, found_intersection)
end

function detect_intersection_inner(grid::HitGrid{N, O}, start::SVector{N, Float64}, dest::SVector{N, Float64}, prev_index::Int32, prev_inside::Bool, voxel::SVector{3, Int}, args...) where {N, O}
    intersection_index = zero(Int32)
    intersection = empty_obstruction_intersections[N]
    for packed in grid.indices[voxel...]
        if grid isa HitGridNoRepeat
            (index, obstruction) = packed
            start_shift = start
            dest_shift = dest
        else
            (index, shift_index, obstruction) = packed
            if iszero(shift_index)
                start_shift = start
                dest_shift = dest
            else
                shift = grid.shifts[shift_index]
                start_shift = start .- shift
                dest_shift = dest .- shift
            end
        end
        if prev_index == index
            new_intersection = detect_intersection(obstruction, start_shift, dest_shift, args..., prev_inside)
        else
            new_intersection = detect_intersection(obstruction, start_shift, dest_shift, args...)
        end
        if (new_intersection.distance >= 0) && (new_intersection.distance < intersection.distance)
            intersection = new_intersection
            intersection_index = index
        end
    end
    return (intersection_index, intersection)
end


function grid_inside_mesh(grid::HitGrid{3, IndexTriangle}, vertices::AbstractVector)
    triangles = map(o -> o.indices, obstructions(grid))
    mean_triangles = map(t->mean(vertices[t]), triangles)
    if grid isa HitGridRepeat
        sz = upper(grid) .- lower(grid)
        mean_triangles = [@. mod(t + sz/2, sz) - sz/2 for t in mean_triangles]
    end
    tree = KDTree(mean_triangles)
    inside_arr = zeros(Bool, size(grid.indices))
    for index in Tuple.(eachindex(IndexCartesian(), inside_arr))
        centre = (@. (index - 0.5) / grid.inv_resolution) .+ lower(grid)

        triangle_index = Int32(nn(tree, centre)[1])
        new_index = triangle_index
        ntry = 0
        while ~iszero(new_index)
            triangle_index = new_index
            (new_index, _) = detect_intersection_grid(grid, mean_triangles[triangle_index], centre, triangle_index, true, vertices)
            ntry += 1
            if ntry > 1000
                break
                error("Grid voxel centre $centre falls exactly on the edge between triangles. Please shift the mesh a tiny amount to fix this.")
            end
        end
        inpr = (centre - mean_triangles[triangle_index]) ⋅ normal(FullTriangle(triangles[triangle_index], vertices))
        if iszero(inpr)
            @warn "Grid voxel centre falls exactly on the mesh element. This will lead to erroneous estimations of what is the inside of a grid. You might want to shift the mesh a tiny amount."
        end
        inside_arr[index...] = inpr < 0
    end
    return inside_arr
end

function isinside(grid::HitGrid{N, IndexTriangle}, position::SVector{N, Float64}, stuck_to::Int32, inside::Bool, vertices, inside_mask) where {N}
    grid_index = get_coordinates(grid, position)
    centre = @. (grid_index - 0.5) / grid.inv_resolution + mesh.grid.lower
    if any(grid_index .< 1) || any(grid_index .> size(mesh.grid.indices))
        return Int32[]
    end

    nhit = 0
    for (index, shift_index, obstruction) in mesh.grid.indices[grid_index...]
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