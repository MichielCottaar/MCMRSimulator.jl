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
BoundingBox(radius::Number) = SymmetricBoundingBox{3}(radius)

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
        PosVector([l[1], l[2], u[3]]),
        PosVector([l[1], u[2], u[3]]),
        PosVector([l[1], u[2], l[3]]),
        PosVector([u[1], l[2], u[3]]),
        PosVector([u[1], u[2], l[3]]),
        PosVector([u[1], l[2], l[3]]),
        u,
    ]
end

function corners(bb::BoundingBox{2}) 
    l = lower(bb)
    u = upper(bb)
    return [
        l,
        SVector{2}([l[1], u[2]]),
        SVector{2}([u[1], l[2]]),
        u,
    ]
end

corners(bb::BoundingBox{1}) = [
    lower(bb),
    upper(bb),
]

isinside(bb::GenericBoundingBox{N}, pos::SVector{N, Float}) where {N} = all(pos .>= lower(bb)) && all(pos .<= upper(bb))
isinside(bb::SymmetricBoundingBox{N}, pos::SVector{N, Float}) where {N} = (maximum(pos) <= bb.radius) && (minimum(pos) >= -bb.radius)


function _bb_intersect_check_dim(::Val{M}, lower::Float, upper::Float, start::SVector{N, Float}, dest::SVector{N, Float}) where {N, M}
    return !(
        (start[M] > upper && dest[M] > upper) ||
        (start[M] < lower && dest[M] < lower) 
    )
end

function _bb_intersect_check_dim(::Val{M}, lower::SVector{N, Float}, upper::SVector{N, Float}, start::SVector{N, Float}, dest::SVector{N, Float}) where {N, M}
    return !(
        (start[M] > upper[M] && dest[M] > upper[M]) ||
        (start[M] < lower[M] && dest[M] < lower[M]) 
    )
end

function _bb_intersect_check_dims(lower, upper, start::SVector{1, Float}, dest::SVector{1, Float})
    return _bb_intersect_check_dim(Val{1}(), lower, upper, start, dest)
end

function _bb_intersect_check_dims(lower, upper, start::SVector{2, Float}, dest::SVector{2, Float})
    return (
        _bb_intersect_check_dim(Val{1}(), lower, upper, start, dest) &&
        _bb_intersect_check_dim(Val{2}(), lower, upper, start, dest)
    )
end

function _bb_intersect_check_dims(lower, upper, start::SVector{3, Float}, dest::SVector{3, Float})
    return (
        _bb_intersect_check_dim(Val{1}(), lower, upper, start, dest) &&
        _bb_intersect_check_dim(Val{2}(), lower, upper, start, dest) &&
        _bb_intersect_check_dim(Val{3}(), lower, upper, start, dest)
    )
end

function possible_intersection(bb::SymmetricBoundingBox{N}, start::SVector{N, Float}, dest::SVector{N, Float}) where {N}
    negative_radius = -bb.radius
    return _bb_intersect_check_dims(negative_radius, bb.radius, start, dest)
end

function possible_intersection(bb::GenericBoundingBox{N}, start::SVector{N, Float}, dest::SVector{N, Float}) where {N}
    return _bb_intersect_check_dims(bb.lower, bb.upper, start, dest)
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

"""
    rotate(bounding_box, rot_mat)

Rotate the bounding box given a (NxM) rotation matrix.
A new bounding box is returned that contains the old one.
"""
function rotate(bb::BoundingBox{N}, rotation_matrix::SMatrix{M, N}) where {N, M}
    all_corners = [rotation_matrix * c for c in corners(bb)]
    lower = SVector{M}([minimum([c[index] for c in all_corners]) for index in 1:M])
    upper = SVector{M}([maximum([c[index] for c in all_corners]) for index in 1:M])
    return BoundingBox(lower, upper)
end