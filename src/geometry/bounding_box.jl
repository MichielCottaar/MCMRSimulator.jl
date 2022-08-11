"""
    BoundingBox(lower::PosVector, upper::PosVector)
    BoundingBox(obstructions::Obstructions)
    BoundingBox(obstruction::Obstruction)

Creates a bounding box containing the [`Obstruction`](@ref) (or [`Obsructions`](@ref)).
For infinitely repeated objects (using [`Repeated`](@ref)) the bounding box of the central object is returned.

Check whether particles are inside using [`isinside`](@ref).
"""
struct BoundingBox
    lower :: PosVector
    upper :: PosVector
    inverse_size :: PosVector
    function BoundingBox(lower::PosVector, upper::PosVector)
        inverse_size = map((l, u) -> isfinite(l) && isfinite(u) ? 1. / (u-l) : 0., lower, upper)
        @assert all(upper .>= lower)
        new(lower, upper, inverse_size)
    end
end

BoundingBox(lower::AbstractVector, upper::AbstractVector) = BoundingBox(PosVector(lower), PosVector(upper))
function BoundingBox(bounding_boxes::BoundingBox...)
    lowers = map(bb->bb.lower, bounding_boxes)
    uppers = map(bb->bb.upper, bounding_boxes)
    BoundingBox(min.(lowers...), max.(uppers...))
end
BoundingBox(obstructions::Obstructions) = BoundingBox(map(o->BoundingBox(o), obstructions)...)

"""
    corners(bb::BoundingBox)

Returns all 8 corners of the [`BoundingBox`](@ref) in an vector of [`PosVector`](@ref).
"""
corners(bb::BoundingBox) = [
    bb.lower,
    SA[bb.lower[1], bb.lower[2], bb.upper[3]],
    SA[bb.lower[1], bb.upper[2], bb.upper[3]],
    SA[bb.lower[1], bb.upper[2], bb.lower[3]],
    SA[bb.upper[1], bb.lower[2], bb.upper[3]],
    SA[bb.upper[1], bb.upper[2], bb.lower[3]],
    SA[bb.upper[1], bb.lower[2], bb.lower[3]],
    bb.upper,
]

isinside(pos::PosVector, bb::BoundingBox) = all(pos .>= bb.lower) && all(pos .<= bb.upper)

"""
    expand(bb::BoundingBox, ratio=1.)

Expand or shrink the [`BoundingBox`](@ref) so that it grows/shrinks by `ratio` along each dimension.
"""
function expand(bb::BoundingBox, ratio=1.)
    sz = bb.upper - bb.lower
    shift = (sz .* ratio .- 1.) .* 0.5
    BoundingBox(bb.lower - shift, bb.upper + shift)
end