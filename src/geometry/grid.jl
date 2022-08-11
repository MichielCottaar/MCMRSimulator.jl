struct RayGridIntersections
    origin :: PosVector
    destination :: PosVector
    direction :: PosVector
end

"""
    ray_grid_intersections([grid, ]origin, destination)

Computes all voxels crossed by a ray between `origin` and `destination` with a [`GridShape`](@ref) (default infinitely extending 1x1x1 grid).
Both origin and destination are length-3 vectors.
The returned object is an iterator returning a tuple with:
- 3-length vector with the voxel that we are crossing through
- Float with the time the ray entered voxel (0=`origin`, 1=`destination`)
- 3-length vector with position within voxel that the ray entered (i.e., numbers between 0 and 1)
- Float with the time the ray left the voxel (0=`origin`, 1=`destination`)
- 3-length vector with position within voxel that the ray left (i.e., numbers between 0 and 1)
"""
ray_grid_intersections(origin :: PosVector, destination :: PosVector) = RayGridIntersections(origin, destination, destination - origin)

Base.iterate(rgi::RayGridIntersections) = Base.iterate(rgi, (rgi.origin, zero(Float), map(o -> Int(floor(o)), rgi.origin)))
function Base.iterate(rgi::RayGridIntersections, state::Tuple{PosVector, Float, SVector{3, Int}})
    (prev_pos::PosVector, prev_time::Float, current_voxel::SVector{3, Int}) = state
    if prev_time >= 1.
        return nothing
    end
    within_voxel = prev_pos - current_voxel
    all_next_hits = map((d, w) -> (d > 0 ? 1. - w : w) / abs(d), rgi.direction, within_voxel)
    dimension = argmin(all_next_hits)
    time_to_hit = all_next_hits[dimension]
    next_time = prev_time + time_to_hit
    if next_time > 1.
        return (
            (current_voxel, prev_time, prev_pos - current_voxel, one(Float), rgi.destination - current_voxel), 
            (rgi.destination, next_time, current_voxel)
        )
    end
    ddim = Int(sign(rgi.direction[dimension]))
    next_voxel = map((i, v) -> i == dimension ? v + ddim : v, 1:3, current_voxel)
    next_pos = rgi.origin .+ (rgi.direction .* next_time)
    res = (
        (current_voxel, prev_time, prev_pos - current_voxel, next_time, next_pos - current_voxel),
        (next_pos, next_time, next_voxel)
    )
    return res
end

Base.length(rgi::RayGridIntersections) = sum(abs.(Int.(floor.(rgi.destination)) .- Int.(floor.(rgi.origin)))) + 1
Base.eltype(::RayGridIntersections) = Tuple{PosVector, Float, PosVector, Float, PosVector}


"""
    GridShape(bounding_box::BoundingBox, size::SVector{3, Int})
    GridShape(bounding_box::BoundingBox, size::Int)

Defines a 3D grid within the [`BoundingBox`](@ref). 
`size` sets the number of voxels along each dimension.
Setting `size` to a single integer will set that number along each dimension 
(except dimensions of infinite size, which are always just one voxel wide).
"""
struct GridShape
    bounding_box :: BoundingBox
    size :: SVector{3, Int}
    voxel_size :: PosVector
    inverse_voxel_size :: PosVector
    function GridShape(bounding_box::BoundingBox, size)
        if isa(size, Int)
            size = map(is -> iszero(is) ? 1 : size, bounding_box.inverse_size)
        end
        size = SVector{3, Int}(size)
        voxel_size = (bounding_box.upper .- bounding_box.lower) ./ size
        @assert all(isfinite.(voxel_size) .|| isone.(size))
        inverse_voxel_size = 1. ./ voxel_size
        new(bounding_box, size, voxel_size, inverse_voxel_size)
    end
end

BoundingBox(g::GridShape) = g.bounding_box
Base.size(g::GridShape) = tuple(g.size)
isinside(pos::PosVector, g::GridShape) = isinside(pos, BoundingBox(g))
function project(pos::PosVector, g::GridShape)
    bb = BoundingBox(g)
    map((p, l, is) -> iszero(is) ? 1.5 : is * (p - l) + 1, pos, bb.lower, g.inverse_voxel_size)
end


function ray_grid_intersections(grid :: GridShape, origin :: PosVector, destination :: PosVector) 
    ray_grid_intersections(project(origin, grid), project(destination, grid))
end