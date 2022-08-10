"""
    Transformed(obstructions::Obstructions, transform::CoordinateTransformations.Transformation)

Underlying `obstructions` are linearly transformed (e.g., rotated or shifted) using the given `transform`.

"""
struct Transformed{N, O, T <: CoordinateTransformations.Transformation} <: Obstruction
    obstrucations :: Obstructions{N, O}
    transform :: T
    inverse :: T
    function Transformed(obstrucations, transform :: CoordinateTransformations.Transformation)
        if isa(obstrucations, Transformed)
            return Transformed(obstrucations.obstrucations, transform ∘ obstrucations.transform )
        else
            o = isa(obstrucations, Obstruction) ? SVector{1}([obstrucations]) : SVector{length(obstrucations)}(obstrucations)
            return new{length(o), eltype(o), typeof(transform)}(o, transform, inv(transform))
        end
    end
end

function detect_collision(movement :: Movement, transform :: Transformed, origin::PosVector) 
    c = detect_collision(
        Movement(
            transform.inverse(movement.origin),
            transform.inverse(movement.destination),
            movement.timestep,
        ),
        transform.obstrucations,
        project(origin, transform)
    )
    if isnothing(c)
        return c
    end
    Collision(
        c.distance,
        transform.transform(c.normal) .- transform.transform(zero(PosVector)),
    )
end

project(pos::PosVector, trans::Transformed) = trans.inverse(pos)
isinside(pos::PosVector, trans::Transformed) = isinside(project(pos, trans), trans.obstrucations)
