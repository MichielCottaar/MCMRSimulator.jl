"""
    Geometry(obstructions...)

Represents the tissue geometry as one or more [`TransformObstruction`](@ref) objects.
Any [`BaseObstruction`](@ref) passed on will be converted into a [`TransformObstruction`](@ref) before being added to the geometry.
"""
struct Geometry{N, T<:TransformObstruction}
    obstructions::SVector{N, T}
end

function Geometry(obstructions::Obstruction...)
    convert_obstruction(o::TransformObstruction) = o
    convert_obstruction(o::BaseObstruction) = TransformObstruction(o)
    Geometry(SVector{length(obstructions)}(convert_obstruction.(obstructions)))
end
Geometry(obstructions::AbstractVector) = Geometry(obstructions...)
Geometry(g::Geometry) = g


isinside(geom::Geometry, pos::PosVector) = any(isinside.(geom.obstructions, pos))

Base.length(geom::Geometry{N}) where {N} = N
    

"""
    detect_collision(movement, geometry[, previous])

Returns a [`Collision`](@ref) object if the given `movement` crosses any obstructions.
The first collision is always returned.
If no collision is detected, `empty_collision` will be returned
"""
function detect_collision(movement :: Movement, geometry::Geometry, previous::Collision=empty_collision)
    collision = empty_collision
    for o in geometry.obstructions
        c_new = detect_collision(movement, o, previous)
        if c_new.distance < collision.distance
            collision = c_new
        end
    end
    collision
end

detect_collision(movement :: Movement, geometry::Geometry{0}, previous::Collision=empty_collision) = empty_collision
detect_collision(movement :: Movement, geometry::Geometry{1}, previous::Collision=empty_collision) = detect_collision(movement, geometry.obstructions[1], previous)

"""
    produces_off_resonance(geometry)

Whether any obstruction produces an off-resonance field.
The field will be computed using [`off_resonance`]
"""
produces_off_resonance(geom::Geometry) = any(produces_off_resonance.(geom.obstructions))

function off_resonance(geom::Geometry{N}, position::PosVector) where {N}
    total = zero(Float)
    for index in 1:N
        obstruction = geom.obstructions[index]
        if produces_off_resonance(obstruction)
            total += off_resonance(obstruction,  position)
        end
    end
    total
end