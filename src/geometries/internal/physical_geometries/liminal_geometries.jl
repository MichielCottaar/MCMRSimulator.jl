"""Physical geometry for a collection of independently placed cell templates."""
module LiminalGeometries

import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, child_type, has_inside, has_single_inside,
    get_intersection_params_requires_inside, inside_indices_eltype, intersection_type,
    inside_indices, find_intersection, find_intersection_requires_inside, get_child
import ..Groups: GeometryTuple
import ..Transformations: Shift
import ...InsideViews: InsideView, child_view

export FixedLiminalGeometry

struct FixedLiminalGeometry{P <: PhysicalGeometry{3}} <: PhysicalGeometry{3}
    geometries::Vector{P}
    number_fractions::Vector{Float64}
    extracellular_fraction::Float64

    function FixedLiminalGeometry(
        geometries::AbstractVector{<:PhysicalGeometry{3}},
        number_fractions::AbstractVector{<:Real},
        extracellular_fraction::Real,
    )
        length(geometries) == length(number_fractions) ||
            throw(ArgumentError("geometries and number_fractions must have the same length"))
        isempty(geometries) && throw(ArgumentError("at least one liminal geometry is required"))
        fractions = Float64[number_fractions...]
        all(isfinite, fractions) && all(>(0), fractions) ||
            throw(ArgumentError("cell number fractions must be finite and positive"))
        fractions ./= sum(fractions)
        child_types = unique(typeof.(geometries))
        child_type = length(child_types) == 1 ?
            only(child_types) : Core.apply_type(Union, child_types...)
        new{child_type}(
            convert(Vector{child_type}, collect(geometries)),
            fractions,
            Float64(extracellular_fraction),
        )
    end
end

function inside_indices(
    ::FixedLiminalGeometry,
    ::SVector{3, Float64},
    intersection=nothing,
)
    throw(ArgumentError(
        "inside_indices is not defined for FixedLiminalGeometry without cell context; " *
        "use the cached inside indices instead",
    ))
end

child_type(::Type{<:FixedLiminalGeometry{P}}) where {P} = Shift{3, P}
has_single_inside(::Type{<:FixedLiminalGeometry}) = false

_child_types(::Type{P}) where {P} = P isa Union ? Base.uniontypes(P) : (P,)
_shifted_child_type(::Type{P}) where {P} = Shift{3, P}
_geometry_tuple_type(::Type{P}) where {P} =
    Core.apply_type(Tuple, (_shifted_child_type(child) for child in _child_types(P))...)

find_intersection_requires_inside(::Type{<:FixedLiminalGeometry}) = Val(true)

function _prepend_liminal_index(::Type{T}) where {T}
    T === Union{} && return Union{}
    T isa Union && return Union{(_prepend_liminal_index(type) for type in Base.uniontypes(T))...}
    T <: Tuple || throw(MethodError(_prepend_liminal_index, (Type{T},)))
    Tuple{Int, SVector{3, Float64}, T.parameters...}
end

inside_indices_eltype(::Type{<:FixedLiminalGeometry{P}}) where {P} =
    _prepend_liminal_index(inside_indices_eltype(GeometryTuple{3, _geometry_tuple_type(P)}))

intersection_type(::Type{<:FixedLiminalGeometry{P}}) where {P} =
    _prepend_liminal_index(intersection_type(GeometryTuple{3, _geometry_tuple_type(P)}))

function get_child(geometry::FixedLiminalGeometry, indices::Tuple)
    index, offset = indices[1:2]
    Shift(geometry.geometries[index], offset), indices[3:end]
end

function find_intersection(
    geometry::FixedLiminalGeometry,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    previous_hit=nothing,
    inside=nothing,
)
    if !isnothing(inside)
        isempty(inside) && throw(ArgumentError("cached inside indices are empty for FixedLiminalGeometry"))
        cached_index = first(inside)
        index, offset = cached_index[1:2]
        child_inside = child_view(inside, (index, offset))
        child_previous = nothing
    elseif !isnothing(previous_hit)
        index, offset = previous_hit[1:2]
        child_inside = nothing
        child_previous = previous_hit[3:end]
    else
        throw(ArgumentError(
            "find_intersection for FixedLiminalGeometry requires cached inside indices " *
            "or a previous intersection",
        ))
    end

    child = Shift(geometry.geometries[index], offset)
    intersection = find_intersection(child, start, destination, child_previous, child_inside)
    isnothing(intersection) ? nothing : (index, offset, intersection...)
end

for trait in (
    :has_inside,
    :get_intersection_params_requires_inside,
)
    @eval function $trait(::Type{<:FixedLiminalGeometry{P}}) where {P}
        $trait(GeometryTuple{3, _geometry_tuple_type(P)})
    end
end

end
