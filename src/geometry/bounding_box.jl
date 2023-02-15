"""
    BoundingBox(lower::PosVector, upper::PosVector)
    BoundingBox(radius::Number)
    BoundingBox(obstruction::Obstruction)
    BoundingBox(obstructions)

Creates a bounding box containing one or more [`Obstruction`](@ref).
For infinitely repeated objects (using [`TransformObstruction`](@ref)) the bounding box of the central object is returned.

Check whether particles are inside using [`isinside`](@ref).

For symmetric obstructions centered at origin (e.g., base cylinders, spheres, walls, annuli, and spirals),
a `SymmetricBoundingBox` is used, which only stores the radius (0 for a wall).
For all other cases a `GenericBoundingBox` is used, which stores the lower and upper bounds as vectors.
"""
abstract type BoundingBox{N} end
 
struct SymmetricBoundingBox{N} <: BoundingBox{N}
    radius :: Float
end

BoundingBox{N}(radius::Number) where {N} = SymmetricBoundingBox{N}(radius)

lower(bb::SymmetricBoundingBox{N}) where {N} = SVector{N, Float}(fill(-bb.radius, N))
upper(bb::SymmetricBoundingBox{N}) where {N} = SVector{N, Float}(fill(bb.radius, N))

struct GenericBoundingBox{N} <: BoundingBox{N}
    lower :: SVector{N, Float}
    upper :: SVector{N, Float}
    inverse_size :: SVector{N, Float}
    function GenericBoundingBox(lower::SVector{N, Float}, upper::SVector{N, Float}) where{N}
        inverse_size = map((l, u) -> isfinite(l) && isfinite(u) ? 1. / (u-l) : 0., lower, upper)
        @assert all(upper .>= lower)
        new{N}(lower, upper, inverse_size)
    end
end

function BoundingBox(lower::AbstractVector, upper::AbstractVector)
    ndim = length(lower)
    GenericBoundingBox(SVector{ndim, Float}(lower), SVector{ndim, Float}(upper))
end

lower(bb::GenericBoundingBox{N}) where {N} = bb.lower
upper(bb::GenericBoundingBox{N}) where {N} = bb.upper

function BoundingBox(bounding_boxes::BoundingBox...)
    lowers = map(lower, bounding_boxes)
    uppers = map(upper, bounding_boxes)
    BoundingBox(min.(lowers...), max.(uppers...))
end

"""
    corners(bb::BoundingBox)

Returns a vector of all corners of the [`BoundingBox`](@ref).
"""
function corners(bb::BoundingBox{3}) 
    l = lower(bb)
    u = upper(bb)
    return [
        l,
        SA[l[1], l[2], u[3]],
        SA[l[1], u[2], u[3]],
        SA[l[1], u[2], l[3]],
        SA[u[1], l[2], u[3]],
        SA[u[1], u[2], l[3]],
        SA[u[1], l[2], l[3]],
        u,
    ]
end

function corners(bb::BoundingBox{2}) 
    l = lower(bb)
    u = upper(bb)
    return [
        l,
        SA[l[1], u[2]],
        SA[u[1], l[2]],
        u,
    ]
end

corners(bb::BoundingBox{1}) = [
    lower(bb),
    upper(bb),
]

isinside(bb::GenericBoundingBox{N}, pos::SVector{N, Float}) where {N} = all(pos .>= lower(bb)) && all(pos .<= upper(bb))
isinside(bb::SymmetricBoundingBox{N}, pos::SVector{N, Float}) where {N} = (maximum(pos) <= bb.radius) && (minimum(pos) >= -bb.radius)

function possible_intersection(bb::SymmetricBoundingBox{N}, start::SVector{N, Float}, dest::SVector{N, Float}) where {N}
    for dim in 1:N
        if (
            (start[dim] > bb.radius && dest[dim] > bb.radius) ||
            (start[dim] < -bb.radius && dest[dim] < -bb.radius) 
        )
            return false
        end
    end
    return true
end

function possible_intersection(bb::GenericBoundingBox{N}, start::SVector{N, Float}, dest::SVector{N, Float}) where {N}
    for dim in 1:N
        if (
            (start[dim] > bb.upper[dim] && dest[dim] > bb.upper[dim]) ||
            (start[dim] < bb.lower[dim] && dest[dim] < bb.lower[dim]) 
        )
            return false
        end
    end
    return true
end

"""
    expand(bb::BoundingBox, ratio=1.)

Expand or shrink the [`BoundingBox`](@ref) so that it grows/shrinks by `ratio` along each dimension.
"""
function expand(bb::BoundingBox, ratio=1.)
    sz = upper(bb) - lower(bb)
    shift = (sz .* ratio .- 1) ./ 2
    BoundingBox(lower(bb) - shift, upper(bb) + shift)
end