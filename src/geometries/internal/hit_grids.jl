module HitGrids

import StaticArrays: SVector
import ..BoundingBoxes: BoundingBox, lower, upper
import ..Obstructions: FixedObstruction, detect_intersection, ObstructionIntersection
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
    shift :: SVector{N, Float64}
end

function (::Type{HitGrid})(obstructions::Vector{<:FixedObstruction{N}}, grid_resolution::Float64, repeats::Union{Nothing, SVector{N, Float64}}=nothing) where {N}
    bounding_boxes = BoundingBox.(obstructions)
    bb_actual = BoundingBox(bounding_boxes)
    if ~isnothing(repeats)
        new_lower = map(bb_actual.lower, bb.actual.upper, repeats) do l, u, r
            if u > r/2
                return -r/2
            else
                return l
            end
        end
        new_upper = map(bb_actual.lower, bb.actual.upper, repeats) do l, u, r
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

    grid = [Tuple{Int32, Int32}[] for _ in Iterators.product(UnitRange.(1, dims)...)]
    for this_bb in bounding_boxes
        if isnothing(repeats)
            range_repeats = SVector{N, UnitRange{Int64}}(fill(N, 0:0))
        else
            range_repeats = UnitRange.(Int.(div.(this_bb.lower, repeats, RoundNearest)), Int.(div.(this_bb.upper, repeats, RoundNearest)))
        end
        for int_shifts in Iterators.product(range_repeats...)
            shift = -SVector{N, Float64}(int_shifts .* repeats)

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
            map(grid_indices) do index_arr
                map(i -> (i[1], obstructions[i[1]]), index_arr)
            end,
        )
    else
        return HitGridRepeat(
            bb,
            dims ./ sz,
            map(grid_indices) do index_arr
                map(i -> (i[1], i[2], obstructions[i[1]]), index_arr)
            end,
            shifts,
        )
    end

end


"""
    get_indices(grid, position)

Get the indices of obstructions that contain the `position`.
"""
function get_indices(grid::HitGrid, position::AbstractVector)
    index = @. Int(div(position - grid.lower, grid.resolution, RoundDown)) + 1
    if any(index .< 1) || any(index .> size(grid.indices))
        return Int32[]
    else
        return grid.indices[index...]
    end
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
        (intersection_index, intersection) = detect_intersection_inner(g, start, dest, prev_index, prev_inside, voxel .+ 1, args...)
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

function detect_intersection_inner(grid::HitGrid{N, O}, start::SVector{N, Float64}, dest::SVector{N, Float64}, prev_index::Int32, prev_inside::Bool, voxel::SVector{3, Int32}, args...) where {N, O}
    intersection_index = zero(Int32)
    intersection = empty_obstruction_intersections[N]
    for packed in grid.indices[voxel...]
        if grid isa HitGridNoRepeat
            (index, obstruction) = packed
            start_shift = start
            dest_shift = dest
        else
            (index, shift, obstruction) = packed
            if iszero(shift)
                start_shift = start
                dest_shift = dest
            else
                shift = g.grid.shifts[shift_index]
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

end