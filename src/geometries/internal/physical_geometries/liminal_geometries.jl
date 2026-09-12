"""Physical geometry for a collection of independently placed cell templates."""
module LiminalGeometries

import StaticArrays: SVector
import Random: rand
import Distributions: Poisson
import ..PhysicalGeometries: PhysicalGeometry, child_type, has_inside, has_single_inside,
    get_intersection_params_requires_inside, inside_indices_eltype, intersection_type,
    inside_indices, find_intersection, find_intersection_requires_inside, get_child,
    volume_sampling, random_surface_positions, bound_intersection_type,
    _merge_types, InternalBoundingBox
import ..Groups: GeometryTuple
import ..Groups: inside_indices_for_any_type
import ..Transformations: Shift
import ...InsideViews: InsideView, child_view
import ...InternalBoundingBoxes
import ...Properties: GeometryProperties

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
)
    throw(ArgumentError(
        "inside_indices is not defined for FixedLiminalGeometry without cell context; " *
        "use the cached inside indices instead",
    ))
end

function inside_indices(
    geometry::FixedLiminalGeometry,
    position::SVector{3, Float64},
    intersection,
)
    cell_index, offset = intersection[1:2]
    child = Shift(geometry.geometries[cell_index], offset)
    child_indices = inside_indices_for_any_type(child, position - offset, intersection[3:end])
    [(cell_index, offset, index...) for index in child_indices]
end

child_type(::Type{<:FixedLiminalGeometry{P}}) where {P} = Shift{3, P}
has_single_inside(::Type{<:FixedLiminalGeometry}) = false

_child_types(::Type{P}) where {P} = P isa Union ? Base.uniontypes(P) : (P,)
_shifted_child_type(::Type{P}) where {P} = Shift{3, P}
_geometry_tuple_type(::Type{P}) where {P} =
    Core.apply_type(GeometryTuple, 3, Core.apply_type(Tuple, (_shifted_child_type(child) for child in _child_types(P))...))
find_intersection_requires_inside(::Type{<:FixedLiminalGeometry}) = Val(true)

function _prepend_liminal_index(::Type{T}) where {T}
    T === Union{} && return Union{}
    T isa Union && return Union{(_prepend_liminal_index(type) for type in Base.uniontypes(T))...}
    T <: Tuple || throw(MethodError(_prepend_liminal_index, (Type{T},)))
    Tuple{Int, SVector{3, Float64}, T.parameters...}
end

inside_indices_eltype(::Type{<:FixedLiminalGeometry{P}}) where {P} =
    _prepend_liminal_index(inside_indices_eltype(_shifted_child_type(P)))

intersection_type(::Type{<:FixedLiminalGeometry{P}}) where {P} =
    _prepend_liminal_index(intersection_type(_shifted_child_type(P)))

function bound_intersection_type(geometry::FixedLiminalGeometry, density::GeometryProperties)
    _merge_types(
        _prepend_liminal_index(bound_intersection_type(child, density.properties[index]))
        for (index, child) in enumerate(geometry.geometries)
    )
end

function get_child(geometry::FixedLiminalGeometry, indices::Tuple)
    index, offset = indices[1:2]
    Shift(geometry.geometries[index], offset), indices[3:end]
end

function _wrapped_offset(
    position::SVector{3, Float64},
    bounding_box::InternalBoundingBox{3},
)
    lower = InternalBoundingBoxes.lower(bounding_box)
    size = 2 .* InternalBoundingBoxes.half_size(bounding_box)
    translated = lower + rand(3) .* size + position
    wrapped = mod.(translated - lower, size) + lower
    wrapped - position
end

function volume_sampling(
    geometry::FixedLiminalGeometry,
    bounding_box::InternalBoundingBox{3},
    volume_density::Number,
)
    volume_density >= 0 || throw(ArgumentError("volume density must be non-negative"))
    intracellular_density = volume_density * (1 - geometry.extracellular_fraction)
    positions = SVector{3, Float64}[]
    indices = Vector{Vector{inside_indices_eltype(typeof(geometry))}}()

    for (cell_index, child) in enumerate(geometry.geometries)
        child_box = InternalBoundingBox(child)
        child_positions, child_indices = volume_sampling(
            child,
            child_box,
            intracellular_density * geometry.number_fractions[cell_index],
        )
        for (position, child_index) in zip(child_positions, child_indices)
            isempty(child_index) && continue
            offset = _wrapped_offset(position, bounding_box)
            push!(positions, position + offset)
            push!(indices, [(cell_index, offset, index...) for index in child_index])
        end
    end

    nsamples = rand(Poisson(geometry.extracellular_fraction * volume_density *
        prod(2 .* InternalBoundingBoxes.half_size(bounding_box))))
    lower = InternalBoundingBoxes.lower(bounding_box)
    size = 2 .* InternalBoundingBoxes.half_size(bounding_box)
    append!(positions, [SVector{3, Float64}(rand(3) .* size .+ lower) for _ in 1:nsamples])
    append!(indices, [Vector{inside_indices_eltype(typeof(geometry))}() for _ in 1:nsamples])
    positions, indices
end

function random_surface_positions(
    geometry::FixedLiminalGeometry,
    density::GeometryProperties,
    bounding_box::InternalBoundingBox{3},
    scale_density,
)
    positions = SVector{3, Float64}[]
    indices = Tuple[]
    intracellular_scale = (1 - geometry.extracellular_fraction)
    for (cell_index, child) in enumerate(geometry.geometries)
        child_positions, child_indices = random_surface_positions(
            child,
            density.properties[cell_index],
            InternalBoundingBox(child),
            scale_density * intracellular_scale * geometry.number_fractions[cell_index],
        )
        for (position, child_index) in zip(child_positions, child_indices)
            offset = _wrapped_offset(position, bounding_box)
            push!(positions, position + offset)
            push!(indices, (cell_index, offset, child_index...))
        end
    end
    positions, indices
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
        $trait(_geometry_tuple_type(P))
    end
end

end
