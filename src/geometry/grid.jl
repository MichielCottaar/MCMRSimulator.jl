struct RayGridIntersections{N}
    origin :: SVector{N, Float}
    destination :: SVector{N, Float}
    direction :: SVector{N, Float}
    abs_inv_direction :: SVector{N, Float}
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
function ray_grid_intersections(origin :: SVector{N, Float}, destination :: SVector{N, Float}) where {N}
    direction = destination - origin
    RayGridIntersections(origin, destination, direction, map(d -> Float(1/abs(d)), direction))
end

Base.iterate(rgi::RayGridIntersections) = Base.iterate(rgi, (rgi.origin, zero(Float), map(o -> Int(floor(o)), rgi.origin)))
function Base.iterate(rgi::RayGridIntersections{N}, state::Tuple{SVector{N, Float}, Float, SVector{N, Int}}) where {N}
    (prev_pos::SVector{N, Float}, prev_time::Float, current_voxel::SVector{N, Int}) = state
    if prev_time >= 1.
        return nothing
    end
    time_to_hit::Float, dimension::Int = next_hit(prev_pos, current_voxel, rgi.direction, rgi.abs_inv_direction)
    next_time = prev_time + time_to_hit
    if next_time > 1.
        return (
            (current_voxel, prev_time, prev_pos - current_voxel, one(Float), rgi.destination - current_voxel), 
            (rgi.destination, next_time, current_voxel)
        )
    end
    ddim = rgi.direction[dimension] > 0 ? 1 : -1
    next_voxel = map((i, v) -> i == dimension ? v + ddim : v, 1:N, current_voxel)
    next_pos = rgi.origin .+ (rgi.direction .* next_time)
    res = (
        (current_voxel, prev_time, prev_pos - current_voxel, next_time, next_pos - current_voxel),
        (next_pos, next_time, next_voxel)
    )
    return res
end

function next_hit(prev_pos::SVector{N, Float}, current_voxel::SVector{N, Int}, direction::SVector{N, Float}, abs_inv_direction::SVector{N, Float}) where {N}
    time_to_hit = Float(Inf)
    dimension = 0
    for dim in 1:N
        within_voxel = prev_pos[dim] - current_voxel[dim]
        d = direction[dim]
        if iszero(d)
            continue
        end
        next_hit = (d > 0. ? 1. - within_voxel : within_voxel) * abs_inv_direction[dim]
        if next_hit < time_to_hit
            time_to_hit = next_hit
            dimension = dim
        end
    end
    return time_to_hit, dimension
end

Base.length(rgi::RayGridIntersections) = sum(abs.(Int.(floor.(rgi.destination)) .- Int.(floor.(rgi.origin)))) + 1
Base.eltype(::RayGridIntersections{N}) where {N} = Tuple{SVector{N, Float}, Float, SVector{N, Float}, Float, SVector{N, Float}}


"""
    GridShape(bounding_box::BoundingBox, size::SVector{3, Int})
    GridShape(bounding_box::BoundingBox, size::Int)

Defines a 3D grid within the [`BoundingBox`](@ref). 
`size` sets the number of voxels along each dimension.
Setting `size` to a single integer will set that number along each dimension 
(except dimensions of infinite size, which are always just one voxel wide).
"""
struct GridShape{N}
    bounding_box :: BoundingBox{N}
    size :: SVector{N, Int}
    voxel_size :: SVector{N, Float}
    inverse_voxel_size :: SVector{N, Float}
    function GridShape(bounding_box::BoundingBox{N}, size) where {N}
        if isa(size, Int)
            size = map(is -> iszero(is) ? 1 : size, bounding_box.inverse_size)
        end
        size = SVector{N, Int}(size)
        voxel_size = (bounding_box.upper .- bounding_box.lower) ./ size
        @assert all(isfinite.(voxel_size) .|| isone.(size))
        inverse_voxel_size = one(Float) ./ voxel_size
        new{N}(bounding_box, size, voxel_size, inverse_voxel_size)
    end
end

BoundingBox(g::GridShape) = g.bounding_box
Base.size(g::GridShape) = tuple(g.size)
isinside(pos::SVector{N, Float}, g::GridShape{N}) where {N} = isinside(pos, BoundingBox(g))
function project(pos::SVector{N, Float}, g::GridShape{N}) where {N}
    bb = BoundingBox(g)
    map((p, l, is) -> iszero(is) ? Float(1.5) : is * (p - l) + 1, pos, bb.lower, g.inverse_voxel_size)
end


function ray_grid_intersections(grid :: GridShape{N}, origin :: SVector{N, Float}, destination :: SVector{N, Float}) where {N}
    ray_grid_intersections(project(origin, grid), project(destination, grid))
end