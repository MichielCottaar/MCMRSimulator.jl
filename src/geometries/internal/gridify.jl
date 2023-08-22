module Gridify

import StaticArrays: SVector
import ..BoundingBoxes: BoundingBox, lower, upper

"""
    Grid(obstructions, vertices, grid_resolution[, repeats])

Predetermines the interactions between the obstructions and the grid.
This is determined purely based on the bounding boxes (could be improved by considering the shape of each obstruction).

Attributes:
- `lower`: lower left coordinate of the final grid (in um).
- `size`: size of the grid (in um).
- `resolution`: resolution of each element in the grid.
- `indices`: Vector of pairs of indices for each element in the grid. The first index indicates the index of the obstruction intersecting with the grid. The second index indicates the shift that should be applied to that obstruction (zero if there is no shift).
- `shifts`: Vector of the shifts that should be applied to obstructions to create a valid repeating grid.
- `repeating`: true if the grid is repeating.
"""
struct Grid{N}
    lower :: SVector{N, Float64}
    size :: SVector{N, Float64}
    resolution :: SVector{N, Float64}
    indices :: Array{Vector{Tuple{Int32, Int32}}, N}
    shifts :: Vector{SVector{N, Float64}}
    repeating :: Bool
end


"""
    get_indices(grid, position)

Get the indices of obstructions that contain the `position`.
"""
function get_indices(grid::Grid, position::AbstractVector)
    index = @. Int(div(position - grid.lower, grid.resolution, RoundDown)) + 1
    if any(index .< 1) || any(index .> size(grid.indices))
        return Int32[]
    else
        return grid.indices[index...]
    end
end


function Grid(obstructions::Vector{BoundingBox{N}}, grid_resolution::Float64, _::Nothing=nothing) where {N}
    bb_actual = BoundingBox(obstructions...)
    extend_by = isfinite(grid_resolution) ? grid_resolution / 100 : 0.001
    bb = BoundingBox(bb_actual.lower .- extend_by, bb_actual.upper .+ extend_by)
    sz = upper(bb) .- lower(bb)
    if isinf(grid_resolution)
        dims = fill(1, N)
    else
        dims = Int.(div.(sz, grid_resolution, RoundUp))
    end
    actual_grid_resolution = sz ./ dims

    grid = [Tuple{Int32, Int32}[] for _ in Iterators.product(UnitRange.(1, dims)...)]
    for (index, this_bb) in enumerate(obstructions)
        l = max.(Int.(div.((lower(this_bb) .- lower(bb)), actual_grid_resolution, RoundUp)), 1)
        u = min.(Int.(div.((upper(this_bb) .- lower(bb)), actual_grid_resolution, RoundDown)) .+ 1, dims)
        for grid_indices in Iterators.product(UnitRange.(l, u)...)
            push!(grid[grid_indices...], (index, 0))
        end
    end
    return Grid{N}(bb.lower, sz, actual_grid_resolution, grid, SVector{N, Float64}[], false)
end

function Grid(obstructions::Vector{BoundingBox{N}}, grid_resolution::Float64, repeats::AbstractVector) where {N}
    repeats = SVector{N}(repeats)
    half_repeats = repeats / 2
    bb = BoundingBox(-half_repeats, half_repeats)

    if isinf(grid_resolution)
        dims = SVector{N, Int}(fill(1, N))
    else
        dims = Int.(div.(repeats, grid_resolution, RoundUp))
    end
    actual_grid_resolution = repeats ./ dims

    grid = [Tuple{Int32, Int32}[] for _ in Iterators.product(UnitRange.(1, dims)...)]
    shifts = SVector{N, Float64}[]
    for (index, this_bb) in enumerate(obstructions)
        range_repeats = UnitRange.(Int.(div.(this_bb.lower, repeats, RoundNearest)), Int.(div.(this_bb.upper, repeats, RoundNearest)))
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
    return Grid(-repeats/2, repeats, actual_grid_resolution, grid, shifts, true)
end

end