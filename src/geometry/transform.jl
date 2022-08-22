"""
    Transformed(obstructions::Obstructions, transform::CoordinateTransformations.Transformation)

Underlying `obstructions` are linearly transformed (e.g., rotated or shifted) using the given `transform`.

"""
struct Transformed{N, O, T <: CoordinateTransformations.Transformation} <: Obstruction
    obstructions :: Obstructions{N, O}
    transform :: T
    inverse :: T
    function Transformed(obstructions, transform :: CoordinateTransformations.Transformation)
        if isa(obstructions, Transformed)
            return Transformed(obstructions.obstructions, transform ∘ obstructions.transform )
        else
            o = isa(obstructions, Obstruction) ? SVector{1}([obstructions]) : SVector{length(obstructions)}(obstructions)
            return new{length(o), eltype(o), typeof(transform)}(o, transform, inv(transform))
        end
    end
end

function detect_collision(movement :: Movement, transform :: Transformed, previous=empty_collision) 
    c = detect_collision(
        Movement(
            transform.inverse(movement.origin),
            transform.inverse(movement.destination),
            movement.timestep,
        ),
        transform.obstructions,
        previous
    )
    if c === empty_collision
        return c
    end
    Collision(
        c.distance,
        transform.transform(c.normal) .- transform.transform(zero(PosVector)),
        c.id,
        c.index
    )
end

project(pos::PosVector, trans::Transformed) = trans.inverse(pos)
isinside(pos::PosVector, trans::Transformed) = isinside(project(pos, trans), trans.obstructions)
function BoundingBox(trans::Transformed)
    input_corners = corners(BoundingBox(trans.obstructions))
    transformed_corners = map(trans.transform, input_corners)
    BoundingBox(min.(transformed_corners...), max.(transformed_corners...))
end

function off_resonance(trans::Transformed, position::PosVector, b0_field::PosVector)
    B0 = project(b0_field, trans) - project(zero(PosVector), trans)
    off_resonance(trans.obstructions, project(position, trans), B0 / norm(B0))
end
