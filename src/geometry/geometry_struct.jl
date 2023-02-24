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


isinside(geom::Geometry, pos::PosVector, stuck_to::Collision) = sum([isinside(o, pos, stuck_to) for o in geom.obstructions])

Base.length(geom::Geometry{N}) where {N} = N
    
function Base.show(io::IO, geom::Geometry{N}) where {N}
    print(io, "Geometry(")
    for obstruction in geom.obstructions
        print(io, obstruction)
        print(io, ", ")
    end
    print(io, ")")
end

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

"""
    inside_MRI_properties(geometry, spin/position, global_props)

Computes the MRI parameters for the spin based on some global settings (`global_props`) and any overriding of those settings in the `geometry`.
"""
function inside_MRI_properties(geom::Geometry{N}, position::PosVector, global_props::MRIProperties) where {N}
    return merge_mri_parameters(SVector{N, MRIProperties}(inside_MRI_properties(o, position) for o in geom.obstructions), global_props)
end
inside_MRI_properties(geom::Geometry{0}, position::PosVector, global_props::MRIProperties) = global_props

inside_MRI_properties(geom::Geometry, position::AbstractVector, global_props::MRIProperties) = inside_MRI_properties(geom, PosVector(position), global_props)

size_scale(geom::Geometry) = minimum(size_scale.(geom.obstructions))
size_scale(geom::Geometry{0}) = Inf

off_resonance_gradient(geom::Geometry) = maximum(off_resonance_gradient.(geom.obstructions))
off_resonance_gradient(geom::Geometry{0}) = zero(Float)

max_timestep_sticking(container::Geometry, default_properties::GlobalProperties, diffusivity::Number) = maximum([max_timestep_sticking(o, default_properties, diffusivity) for o in container.obstructions])
max_timestep_sticking(container::Geometry{0}, default_properties::GlobalProperties, diffusivity::Number) = Inf

function random_surface_spins(geometry::Geometry, bounding_box::BoundingBox, volume_density::Number, default_surface_density::Number; kwargs...)
    vcat([random_surface_spins(o, bounding_box, volume_density, default_surface_density; kwargs...) for o in geometry.obstructions]...)
end