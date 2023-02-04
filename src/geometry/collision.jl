"""
    Collision(distance, normal, properties; index=0, inside=false)

A detected collision along the movement.

# Parameters
- `distance`: number between 0 and 1 indicating the distance of the collision from the origin (0) to the destination point (1)
- `normal`: normal of the obstruction at the collision site. To get correct reflection the normal should point in the direction of the incoming particle.
- `properties`: [`ObstructionProperties`](@ref) of the obstruction the spin collided with.
- `index`: Index of which triangle in [`Mesh`](@ref) got hit
- `inside`: Whether the obstruction was hit on the inside or the outside. For a mesh triangle the outside is considered in the direction of the normal. For a wall the outside is in the positive direction.
"""
struct Collision
    distance :: Float
    normal :: PosVector
    properties :: ObstructionProperties
    index :: Int
    inside :: Bool
    Collision(distance, normal, properties, index, inside) = new(iszero(distance) ? distance : prevfloat(distance), normal, properties, index, inside)
end

Collision(distance, normal, properties; index=0, inside=false) = Collision(distance, normal, properties, index, inside)
ObstructionProperties(c :: Collision) = c.properties


const empty_collision = Collision(Inf, SA[0, 0, 0], ObstructionProperties(), 0, false)

for accessor in (:MT_fraction, :permeability)
    @eval function $(accessor)(c::Collision, defaults) 
        value = $(accessor)(c.properties)
        if isnan(value)
            value = $(accessor)(defaults)
        end
        @assert !isnan(value)
        return value
    end
end

