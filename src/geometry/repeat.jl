"""
    Repeated(obstructions, repeats)

Underlying `obstructions` are repeated ad infinitum as described in `repeats` (a length-3 vector).
Zeros or infinity in `repeats` indicate that no repeating should be applied in this direction.
To repeat in a different direction that the cardinal directions first apply `Repeated` and then [`Transformed`](@ref).

If the obstructions are larger than the repeat size they will be cut off!
"""
struct Repeated{N, T} <: Obstruction
    obstructions :: Obstructions{N, T}
    repeats :: PosVector
    function Repeated(obstructions, repeats)
        o = isa(obstructions, Obstruction) ? SVector{1}([obstructions]) : SVector{length(obstructions)}(obstructions)
        rs = SVector{3, Float}([iszero(r) ? Inf : abs(r) for r in repeats])
        any(isfinite.(rs)) ? new{length(o), eltype(o)}(o, rs) : obstructions
    end
end

project(pos::PosVector, repeat::Repeated) = map((p, r) -> isfinite(r) ? p - div(p, r, RoundNearest) * r : p, pos, repeat.repeats)
function BoundingBox(repeat::Repeated)
    bb = BoundingBox(repeat.obstructions)
    BoundingBox(
        map((l, r) -> isfinite(r) ? max(l, -r/2) : l, bb.lower, repeat.repeats),
        map((u, r) -> isfinite(r) ? min(u, r/2) : u, bb.upper, repeat.repeats),
    )
end

function detect_collision(movement :: Movement, repeat :: Repeated, previous=empty_collision)
    fdiv(p, r) = p / r + 1//2
    origin = map(fdiv, movement.origin, repeat.repeats)
    destination = map(fdiv, movement.destination, repeat.repeats)
    for (_, t1, p1, t2, p2) in ray_grid_intersections(origin, destination)
        f(r, p, p_orig) = isfinite(r) ? r * (p - 1//2) : p_orig
        pos1 = f.(repeat.repeats, p1, (movement.destination .* t1) .+ (movement.origin .* (1 - t1)))
        pos2 = f.(repeat.repeats, p2, (movement.destination .* t2) .+ (movement.origin .* (1 - t2)))
        c = detect_collision(
            Movement(pos1, pos2, Float(1)),
            repeat.obstructions,
            previous
        )
        if c !== empty_collision
            return Collision(
                c.distance * (t2 - t1) + t1,
                c.normal,
                c.id,
                c.index
            )
        end
    end
    return empty_collision
end


isinside(pos::PosVector, repeat::Repeated) = isinside(project(pos, repeat), repeat.obstructions)
