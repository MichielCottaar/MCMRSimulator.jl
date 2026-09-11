"""Physical geometry for a collection of independently placed cell templates."""
module LiminalGeometries

import ..PhysicalGeometries: PhysicalGeometry

export FixedLiminalGeometry

struct FixedLiminalGeometry{P <: PhysicalGeometry{3}} <: PhysicalGeometry{3}
    geometries::Vector{P}
    volume_fractions::Vector{Float64}
    extracellular_fraction::Float64

    function FixedLiminalGeometry(
        geometries::AbstractVector{<:PhysicalGeometry{3}},
        volume_fractions::AbstractVector{<:Real},
        extracellular_fraction::Real,
    )
        length(geometries) == length(volume_fractions) ||
            throw(ArgumentError("geometries and volume_fractions must have the same length"))
        isempty(geometries) && throw(ArgumentError("at least one liminal geometry is required"))
        child_types = unique(typeof.(geometries))
        child_type = length(child_types) == 1 ?
            only(child_types) : Core.apply_type(Union, child_types...)
        new{child_type}(
            convert(Vector{child_type}, collect(geometries)),
            Float64[volume_fractions...],
            Float64(extracellular_fraction),
        )
    end
end

end
