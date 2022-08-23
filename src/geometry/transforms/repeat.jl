"""
    Repeated(obstruction, repeats)

Underlying `obstruction` is repeated ad infinitum as described in `repeats` (a length-3 vector).
Zeros or infinity in `repeats` indicate that no repeating should be applied in this direction.

If the obstruction is larger than the repeat size they will be cut off!
"""
struct Repeated{O<:BaseObstruction{N}, } <: TransformObstruction{N}
    obstruction :: O
    repeats :: SVector{N, Float}
    lorentz_radius :: Float
    lorentz_repeats :: SVector{N, Int}
    chi
    function Repeated(obstruction::BaseObstruction{N}, repeats; lorentz_radius=1.) where {N}
        rs = SVector{N, Float}([iszero(r) ? Inf : abs(r) for r in repeats])
        lorentz_repeats = map(r -> isfinite(r) ? Int(div(lorentz_radius, r, RoundUp)) : 0, rs)
        chi_mult = N == 2 ? 2 * π : 4 * π / 3
        chi = chi_mult * total_susceptibility(obstruction) / prod(rs)
        new{length(o), eltype(o)}(o, rs, lorentz_radius, lorentz_repeats, chi)
    end
end

project(repeat::Repeated{N}, pos::SVector{N, Float}) where {N} = map((p, r) -> isfinite(r) ? p - div(p, r, RoundNearest) * r : p, pos, repeat.repeats)

function BoundingBox(repeat::Repeated)
    bb = BoundingBox(repeat.obstructions)
    BoundingBox(
        map((l, r) -> isfinite(r) ? max(l, -r/2) : l, bb.lower, repeat.repeats),
        map((u, r) -> isfinite(r) ? min(u, r/2) : u, bb.upper, repeat.repeats),
    )
end

function detect_collision(movement :: Movement{N}, repeat :: Repeated{N}, previous=empty_collision) where {N}
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


function off_resonance(obstruction::Repeated{N}, position::SVector{N, Float}, b0_field::PosVector) where {N}
    # Contribution from outside of the Lorentz cavity
    if N == 2
        field = obstruction.chi * (1 - b0_field[3] * b0_field[3])
    else
        field = obstruction.chi
    end
    # Contribution from within the Lorentz cavity
    return field + repeated_off_resonance(obstruction.obstruction, position, b0_field, obstruction.repeats, obstruction.lorentz_radius, obstruction.lorentz_repeats)
end