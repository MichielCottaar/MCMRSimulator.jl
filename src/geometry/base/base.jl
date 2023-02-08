"""Base obstruction type hindering free diffusion.

These obstructions are always perfectly aligned with the main axes and centered at the origin (except for [`Mesh`](@ref)).
They can be moved/rotated/repeated by applying [`TransformObstruction`](@ref).

The dimensionality N indicates the dimensionality of the input data 
(1 for [`Wall`](@ref), 2 for [`Cylinder`](@ref), 3 for [`Sphere`](@ref) or [`Mesh`](@ref)).

Each obstruction needs to define the following interface:
- [`detect_collision`](@ref)(movement, obstruction, previous_collision): returns any interesection between the movement and the obstruction
- [`produces_off_resonance`](@ref)(obstruction), optional: whether the obstruction produces an off-resonance field (false by default). If true, the obstruction should also define:
    - [`lorentz_off_resonance`](@ref)(obstruction, position, ...): computes the off-resonance due to the obstruction at position
    - [`total_susceptibility`](@ref)(obstruction): computes total susceptibility caused by this obstruction
    - [`off_resonance_gradient`](@ref)(obstruction): computes the maximum off-resonance gradient induced by this obstruction
- [`isinside`](@ref)(obstruction, position): true if position is inside the obstruction
- [`BoundingBox`](@ref)(obstruction): returns a [`BoundingBox`] fully containing the obstruction
- [`size_scale`](@ref)(obstruction): returns a scale of the size of the object; used to estimate a maximum size step in the simulation
"""
abstract type BaseObstruction{N} <: Obstruction{N} end

ObstructionProperties(obstruction :: BaseObstruction) = obstruction.properties

"""
    collided(o::BaseObstruction, c::Collision)

Returns true if the collision `c` hit the obstruction `o`.
"""
collided(o::BaseObstruction, c::Collision) = ObstructionProperties(o) == ObstructionProperties(c)

"""
    inside_MRI_properties(obstruction, position)

Returns the inside MRI properties of the obstruction only if the position is within the property.
For positions outside of the MRI properties an empty [`MRIProperties`](@ref) is returned.
"""
function inside_MRI_properties(obstruction::BaseObstruction{N}, position::SVector{N, Float}) where{N}
    if empty_mri_properties(obstruction.properties.inside) || isinside(obstruction, position)
        return obstruction.properties.inside
    else
        return MRIProperties()
    end
end

# Base obstruction API
"""
    detect_collision(movement, base_obstruction, previous_collision)

Returns any intersection between the [`Movement`](@ref) and the [`BaseObstruction`].
"""
function detect_collision end

"""
    produces_off_resonance(base_obstruction)

Whether the obstruction produces an off-resonance field.
The field will be computed using [`lorentz_off_resonance`]
"""
produces_off_resonance(obstruction::BaseObstruction) = false


"""
    lorentz_off_resonance(base_obstruction, position, b0_field, repeat_dist, radius, nrepeats)

Computes the off-resonance field produced by the obstruction ([`BaseObstruction`](@ref)) at the given `position` given the `b0_field`.
The field generated by any repeats of the base obstruction within a distance of `radius` will also be considered (presuming the obstruction repeats every `repeat_dist`).
This maximum number of repeats that need to be considered is precomputed as `nrepeats`.
"""
function lorentz_off_resonance end


"""
    isinside(obstruction/geometry/bounding_box, position)
    isinside(obstructions/geometry/bounding_box, spin)
    isinside(obstructions/geometry/bounding_box, snapshot)

Test whether the particles are inside a [`BaseObstruction`](@ref), [`TransformObstruction`](@ref), [`Geometry`](@ref) or [`BoundingBox`](@ref) object.
"""
isinside(something, pos::AbstractVector) = isinside(something, SVector{length(pos)}(pos))
isinside(something, spin::Spin) = isinside(something, spin.position)
isinside(something, snapshot::Snapshot) = [isinside(something, spin) for spin in snapshot]

"""
    total_susceptibility(obstruction)

Computes total surface (for 2D) or volume (for 3D) susceptibility of a base obstruction.
"""
total_susceptibility(o::BaseObstruction) = zero(Float)

"""
    off_resonance_gradient(obstruction[, B0])
    off_resonance_gradient(geometry[, B0])

Computes the maximum off-resonance gradient produced by this obstruction.
For a geometry with multiple obstructions, the maximum value is returned.
If a B0 field is provided the result is returned in kHz/um, otherwise in ppm/um.
"""
off_resonance_gradient(obstruction::BaseObstruction) = zero(Float)
off_resonance_gradient(obstruction_or_geometry, B0::Number) = off_resonance_gradient(obstruction_or_geometry) * 1e-6 * B0 * gyromagnetic_ratio


"""
    size_scale(geometry)
    size_scale(obstruction)

Computes the minimum size scale of the obstructions in the geometry.

The size scale for each obstruction is defined as:
- [`cylinders`](@ref)/[`spheres`](@ref): minimum radius
- [`annuli`](@ref): minimum inner radius
- [`mesh`](@ref): square root of median triangle size
- [`wall`](@ref): distance between closest walls (infinite for single, non-repeating walls)
"""

include("random_draws.jl")
include("cylinder.jl")
include("sphere.jl")
include("wall.jl")
include("mesh.jl")
include("annulus.jl")
include("spiral.jl")