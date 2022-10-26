"""
Supertype of any obstruction to the free diffusion of water. Some might also generate off-resonance fields.

There are two types of obstructions:
- [`BaseObstruction`](@ref) with the basic obstructions (walls, sphere, cylinders, meshes)
- [`TransformObstruction`](@ref), which transform the base obstructions
"""
abstract type Obstruction{N} end

"""
    isinside(obstructions/bounding_box, position)
    isinside(obstructions/bounding_box, spin)
    isinside(obstructions/bounding_box, snapshot)

Test whether the particles are inside any of the [`Obstruction`](@ref) objects (or in the [`BoundingBox`](@ref)).
"""
isinside(o, pos::Vector) = isinside(o, PosVector(pos))
isinside(o, spin::Spin) = isinside(o, position(spin))
isinside(o, snapshot::Snapshot) = map(s -> isinside(o, s), snapshot.spins)
isinside(obstructions::Tuple, pos::PosVector) = maximum(o -> isinside(o, pos), obstructions)
isinside(obstructions::AbstractVector{<:Obstruction}, position::PosVector) = isinside(Tuple(obstructions), position)


"""
    project(position, transform::TransformObstruction)

Computes the position in the space of the obstructions wrapped by the [`TransformObstruction`](@ref).

    project(position, grid::GridShape)

Computes the voxel index for the position on the [`GridShape`](@ref). 
This will return a result even if the point is outside of the grid. Use [`isinside`](@ref)(position, grid) to check that.
"""
function project end


"""
    off_resonance(obstructions, position[, b0_field])

Computes the off-resonance field at the position due to the obstructions with the magnetic field orientation from `b0_field`.
"""
function off_resonance(obstructions::Tuple, position::PosVector, b0_field=PosVector([0, 0, 1])::PosVector)
    total = zero(Float)
    for o in obstructions
        total += off_resonance(o, position, b0_field)
    end
    total
end
off_resonance(obstructions::Obstruction, position::PosVector, b0_field::PosVector) = zero(Float)
off_resonance(obstructions::Tuple{}, position::PosVector, b0_field=PosVector([0, 0, 1])::PosVector) = zero(Float)
off_resonance(obstructions::AbstractVector{<:Obstruction}, position::PosVector, b0_field=PosVector([0, 0, 1])::PosVector) = sum(o->off_resonance(o, position, b0_field), obstructions)

"""
    lorentz_off_resonance(obstructions, position, b0_field, repeat_dist, radius, nrepeats)

Computes the off-resonance field contribution of repeating compartments within a spherical or cylindrical Lorentz cavity.
"""
lorentz_off_resonance(obstruction::Obstruction{N}, position::SVector{N, Float}, b0_field::SVector{N, Float}, repeat_dist::SVector{N, Float}, radius::Float, nrepeats::SVector{N, Int}) where {N} = zero(Float)


"""
    total_susceptibility(obstruction)

Computes the total magnetic susceptibility of an obstructions.
This might be integrated over the surface or the volume depending on the dimensionality of the obstruction.
"""
total_susceptibility(obstruction::Obstruction) = zero(Float)

include("bounding_box.jl")
include("properties.jl")
include("diffuse.jl")
include("grid.jl")
include("base/base.jl")
include("transform.jl")