"""Physical geometry for a collection of independently placed cell templates."""
module LiminalGeometries

import ..PhysicalGeometries: PhysicalGeometry, child_type, has_inside, has_single_inside,
    get_intersection_params_requires_inside, inside_indices_eltype, intersection_type
import ..Groups: GeometryTuple

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

child_type(::Type{<:FixedLiminalGeometry{P}}) where {P} = P
has_single_inside(::Type{<:FixedLiminalGeometry}) = false

_child_types(::Type{P}) where {P} = P isa Union ? Base.uniontypes(P) : (P,)
_geometry_tuple_type(::Type{P}) where {P} = Core.apply_type(Tuple, _child_types(P)...)

for trait in (
    :has_inside,
    :get_intersection_params_requires_inside,
    :inside_indices_eltype,
    :intersection_type,
)
    @eval function $trait(::Type{<:FixedLiminalGeometry{P}}) where {P}
        $trait(GeometryTuple{3, _geometry_tuple_type(P)})
    end
end

end
