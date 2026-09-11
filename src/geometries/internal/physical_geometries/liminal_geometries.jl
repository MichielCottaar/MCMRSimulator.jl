"""Physical geometry for a collection of independently placed cell templates."""
module LiminalGeometries

import ..PhysicalGeometries: PhysicalGeometry

export FixedLiminalGeometry

struct FixedLiminalGeometry <: PhysicalGeometry{3}
    geometries::Vector{PhysicalGeometry{3}}
    volume_fractions::Vector{Float64}
    extracellular_fraction::Float64

    function FixedLiminalGeometry(
        geometries::AbstractVector{<:PhysicalGeometry{3}},
        volume_fractions::AbstractVector{<:Real},
        extracellular_fraction::Real,
    )
        length(geometries) == length(volume_fractions) ||
            throw(ArgumentError("geometries and volume_fractions must have the same length"))
        new(
            PhysicalGeometry{3}[geometries...],
            Float64[volume_fractions...],
            Float64(extracellular_fraction),
        )
    end
end

end
