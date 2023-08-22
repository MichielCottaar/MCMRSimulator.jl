"""
Types:
- [`Susceptibility`](@ref)

Methods:
- [`single_susceptibility`](@ref)
"""
module Base

"""
    BaseSusceptibility{N}

Parent type of all susceptibility sources.
Each `BaseSusceptibility` object represents a single susceptibility source in `N`-dimensional space 
(e.g., N=2 for [`CylinderSusceptibility`](@ref) or N=3 for [`PointSusceptibility`](@ref)).
"""
abstract type BaseSusceptibility{N} end

"""
    single_susceptibility(source::BaseSusceptibility, offset, distance, stuck_to)

Computes the off-resonance field contribution from a [`BaseSusceptibility`](@ref) source given:
- `offset` vector from spin position to the centre of the `source`
- `distance` of the spin from the source (which is the norm of `offset`)
- `stuck_to` tuple with the pair of indices that the particle is stuck to as well as a boolean indicating whether it is stuck to the inside.
"""
function single_susceptibility end

"""
    single_susceptibility_gradient(source::BaseSusceptibility)

Maximal off-resonance gradient in ppm/um produced by a single susceptibility source.
"""
function single_susceptibility_gradient end


end
