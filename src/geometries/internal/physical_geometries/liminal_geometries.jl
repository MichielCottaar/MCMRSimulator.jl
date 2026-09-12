"""Physical geometry for a collection of independently placed cell templates."""
module LiminalGeometries

import StaticArrays: SVector
import LinearAlgebra: norm, ⋅
import Random: rand, randperm
import Distributions: Poisson
import ..PhysicalGeometries: PhysicalGeometry, child_type, has_inside, has_single_inside,
    get_intersection_params_requires_inside, inside_indices_eltype, intersection_type,
    inside_indices, find_intersection, find_intersection_requires_inside, get_child,
    volume_sampling, random_surface_positions, bound_intersection_type,
    _merge_types, InternalBoundingBox, estimate_surface, estimate_volume,
    get_intersection_params, to_inside_index,
    projected_surface_area, inverse_mean_free_path
import ..Groups: GeometryTuple
import ..Groups: inside_indices_for_any_type, _append_type
import ..Transformations: Shift
import ...InsideViews: InsideView, child_view
import ...InternalBoundingBoxes
import ....BoundingBoxes: BoundingBoxNotSupported
import ...Properties: GeometryProperties, GeometryLeafProperties

export FixedLiminalGeometry, OuterSurfaceSampling, sample!

"""Mutable library of sampled outer-surface data for liminal cells.

`surface_indices` stores raw intersection indices followed by the sampled
inside-side flag. Samples are consumed in order starting at
`index_to_sample`.
"""
mutable struct OuterSurfaceSampling{I}
    positions::Vector{SVector{3, Float64}}
    normals::Vector{SVector{3, Float64}}
    cell_indices::Vector{Int}
    surface_indices::Vector{I}
    weight::Float64
    index_to_sample::Int
    lock::ReentrantLock
end

function projected_surface_area(
    sampling::OuterSurfaceSampling,
    direction::SVector{3, Float64},
)
    direction_norm = norm(direction)
    iszero(direction_norm) && throw(ArgumentError("direction must be non-zero"))
    unit_direction = direction / direction_norm
    projected_sum = sum(
        abs(normal ⋅ unit_direction)
        for normal in sampling.normals
    ) / 2
    sampling.weight * projected_sum
end

struct FixedLiminalGeometry{P <: PhysicalGeometry{3}} <: PhysicalGeometry{3}
    geometries::Vector{P}
    number_fractions::Vector{Float64}
    extracellular_fraction::Float64
    total_surface_area::Float64
    weighted_cell_volume::Float64
    outer_surface_sampling::OuterSurfaceSampling

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
        total_surface_area = sum(
            fraction * estimate_surface(child; outer=true).area
            for (fraction, child) in zip(fractions, geometries)
        )
        weighted_cell_volume = sum(
            fraction * estimate_volume(child)
            for (fraction, child) in zip(fractions, geometries)
        )
        child_types = unique(typeof.(geometries))
        child_type = length(child_types) == 1 ?
            only(child_types) : Core.apply_type(Union, child_types...)
        surface_index_type = _append_type(
            _merge_types(intersection_type(child) for child in _child_types(child_type)),
            Bool,
        )
        fixed = new{child_type}(
            convert(Vector{child_type}, collect(geometries)),
            fractions,
            Float64(extracellular_fraction),
            total_surface_area,
            weighted_cell_volume,
            OuterSurfaceSampling{surface_index_type}(
                SVector{3, Float64}[],
                SVector{3, Float64}[],
                Int[],
                surface_index_type[],
                0.,
                1,
                ReentrantLock(),
            ),
        )
        sample!(fixed.outer_surface_sampling, fixed)
        fixed
    end
end

InternalBoundingBox(::FixedLiminalGeometry) = throw(BoundingBoxNotSupported(
    "liminal geometries do not have a finite bounding box",
))

function inverse_mean_free_path(
    geometry::FixedLiminalGeometry,
    direction::SVector{3, Float64},
)
    iszero(geometry.extracellular_fraction) && return Inf
    cell_number_density = (1 - geometry.extracellular_fraction) /
        (geometry.extracellular_fraction * geometry.weighted_cell_volume)
    cell_number_density * projected_surface_area(geometry.outer_surface_sampling, direction)
end

function _is_outer_surface_sample(geometry::PhysicalGeometry, position, full_index)
    collision_indices = full_index[1:(end - 1)]
    current_inside = to_inside_index(geometry, collision_indices)
    all(
        inside_index == current_inside
        for inside_index in inside_indices_for_any_type(geometry, position, nothing)
    )
end

function sample!(
    sampling::OuterSurfaceSampling,
    geometry::FixedLiminalGeometry,
    N::Integer=1000,
)
    N >= 0 || throw(ArgumentError("number of surface samples must be non-negative"))
    empty!(sampling.positions)
    empty!(sampling.normals)
    empty!(sampling.cell_indices)
    empty!(sampling.surface_indices)
    sampling.weight = 0.
    sampling.index_to_sample = 1
    iszero(N) && return sampling
    iszero(geometry.total_surface_area) && return sampling

    sample_density = N / geometry.total_surface_area
    sampling.weight = inv(sample_density)
    density = GeometryLeafProperties(1.0)
    for (cell_index, child) in enumerate(geometry.geometries)
        scale_density = sample_density * geometry.number_fractions[cell_index]
        positions, indices = random_surface_positions(
            child,
            density,
            InternalBoundingBox(child),
            scale_density,
        )
        for (position, full_index) in zip(positions, indices)
            _is_outer_surface_sample(child, position, full_index) || continue
            params = get_intersection_params(
                child,
                position,
                position,
                (full_index..., 0.0),
            )
            normal = full_index[end] ? -params.normal : params.normal
            push!(sampling.positions, position)
            push!(sampling.normals, normal)
            push!(sampling.cell_indices, cell_index)
            push!(sampling.surface_indices, full_index)
        end
    end
    permutation = randperm(length(sampling.positions))
    sampling.positions .= sampling.positions[permutation]
    sampling.normals .= sampling.normals[permutation]
    sampling.cell_indices .= sampling.cell_indices[permutation]
    sampling.surface_indices .= sampling.surface_indices[permutation]
    sampling
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
    isnothing(intersection) && throw(ArgumentError(
        "inside_indices is not defined for FixedLiminalGeometry without cell context; " *
        "use the cached inside indices instead",
    ))
    cell_index, offset = intersection[1]
    child = Shift(geometry.geometries[cell_index], offset)
    child_indices = inside_indices_for_any_type(child, position - offset, intersection[2:end])
    [tuple((cell_index, offset), index...) for index in child_indices]
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
    Tuple{Tuple{Int, SVector{3, Float64}}, T.parameters...}
end

inside_indices_eltype(::Type{<:FixedLiminalGeometry{P}}) where {P} =
    _prepend_liminal_index(
        _merge_types(inside_indices_eltype(_shifted_child_type(child)) for child in _child_types(P)),
    )

intersection_type(::Type{<:FixedLiminalGeometry{P}}) where {P} =
    _prepend_liminal_index(
        _merge_types(intersection_type(_shifted_child_type(child)) for child in _child_types(P)),
    )

function bound_intersection_type(geometry::FixedLiminalGeometry, density::GeometryProperties)
    _merge_types(
        _prepend_liminal_index(bound_intersection_type(child, density.properties[index]))
        for (index, child) in enumerate(geometry.geometries)
    )
end

function get_child(geometry::FixedLiminalGeometry, indices::Tuple)
    index, offset = indices[1]
    Shift(geometry.geometries[index], offset), indices[2:end]
end

function to_property_index(geometry::FixedLiminalGeometry, indices::Tuple)
    child, child_indices = get_child(geometry, indices)
    cleaned = to_property_index(child, child_indices)
    (indices[1][1], cleaned...)
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
    density_scale = prod(2 .* InternalBoundingBoxes.half_size(bounding_box)) /
        geometry.weighted_cell_volume
    positions = SVector{3, Float64}[]
    indices = Vector{Vector{inside_indices_eltype(typeof(geometry))}}()

    for (cell_index, child) in enumerate(geometry.geometries)
        child_box = InternalBoundingBox(child)
        child_positions, child_indices = volume_sampling(
            child,
            child_box,
            intracellular_density * geometry.number_fractions[cell_index] * density_scale,
        )
        for (position, child_index) in zip(child_positions, child_indices)
            isempty(child_index) && continue
            offset = _wrapped_offset(position, bounding_box)
            push!(positions, position + offset)
            push!(indices, [tuple((cell_index, offset), index...) for index in child_index])
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
            push!(indices, tuple((cell_index, offset), child_index...))
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
    displacement = destination - start
    distance = norm(displacement)
    iszero(distance) && return nothing
    liminal_step = (isnothing(inside) || isempty(inside)) &&
        (isnothing(previous_hit) || !previous_hit[end - 1])

    if liminal_step
        direction = displacement / distance
        inverse_path = inverse_mean_free_path(geometry, direction)
        iszero(inverse_path) && return nothing
        encounter_distance = -log1p(-rand()) / inverse_path
        encounter_distance > distance && return nothing
        sampling = geometry.outer_surface_sampling
        lock(sampling.lock)
        sample_position = zero(SVector{3, Float64})
        sample_cell_index = 0
        sample_surface_index = ()
        try
            selected = 0
            while true
                if sampling.index_to_sample > length(sampling.positions)
                    sample!(sampling, geometry)
                end
                selected = sampling.index_to_sample
                sampling.index_to_sample += 1
                acceptance_probability = max(0., -direction ⋅ sampling.normals[selected])
                if rand() < acceptance_probability
                    break
                end
            end
            sample_position = sampling.positions[selected]
            sample_cell_index = sampling.cell_indices[selected]
            sample_surface_index = sampling.surface_indices[selected]
        finally
            unlock(sampling.lock)
        end
        encounter_position = start + encounter_distance * direction
        offset = encounter_position - sample_position
        return (
            (
                sample_cell_index,
                offset,
            ),
            sample_surface_index[1:(end - 1)]...,
            false,
            encounter_distance / distance,
        )
    elseif !isnothing(inside)
        isempty(inside) && throw(ArgumentError("cached inside indices are empty for FixedLiminalGeometry"))
        cached_index = first(inside)
        index, offset = cached_index[1]
        child_inside = child_view(inside, (cached_index[1],))
        child_previous = isnothing(previous_hit) ? nothing : previous_hit[2:end]
    elseif !isnothing(previous_hit)
        index, offset = previous_hit[1]
        child_inside = nothing
        child_previous = previous_hit[2:end]
    else
        throw(ArgumentError("invalid FixedLiminalGeometry intersection state"))
    end

    child = Shift(geometry.geometries[index], offset)
    intersection = find_intersection(child, start, destination, child_previous, child_inside)
    isnothing(intersection) ? nothing : ((index, offset), intersection...)
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
