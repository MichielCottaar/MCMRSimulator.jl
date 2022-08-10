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

isinside(pos::PosVector, bb::BoundingBox) = all(pos .> bb.lower) && all(pos .< bb.upper)

