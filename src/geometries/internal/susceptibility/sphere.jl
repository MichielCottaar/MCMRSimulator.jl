"""
Off-resonance field due to the constant magnetic susceptibility of a sphere.
"""
module Sphere

import StaticArrays: SVector
import ..Base: BaseSusceptibility, single_susceptibility, single_susceptibility_gradient

"""
    SphereSusceptibility(radius, susceptibility)

Creates a susceptibility source for a sphere with constant scalar susceptibility.
The field inside the sphere is constant; the field outside is represented exactly
by the dipole approximation used by the susceptibility grid.
"""
struct SphereSusceptibility <: BaseSusceptibility{3}
    radius :: Float64
    internal_field :: Float64
    function SphereSusceptibility(radius::Number, susceptibility::Number)
        new(Float64(radius), Float64(susceptibility) / 3)
    end
end

function single_susceptibility(
    sphere::SphereSusceptibility,
    position::AbstractVector,
    distance::Number,
    stuck_inside::Union{Nothing, Bool},
    b0_field::SVector{3, Float64},
)
    sphere.internal_field
end

single_susceptibility_gradient(sphere::SphereSusceptibility) =
    6 * abs(sphere.internal_field) / sphere.radius

end
