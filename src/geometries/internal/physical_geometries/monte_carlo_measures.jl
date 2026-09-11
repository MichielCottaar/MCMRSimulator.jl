module MonteCarloMeasures

import Random: AbstractRNG, default_rng, rand
import StaticArrays: SVector
import ..PhysicalGeometries: PhysicalGeometry, IntersectionParams, find_intersection,
    get_intersection_params, random_surface_positions, has_inside, has_single_inside,
    inside_indices, isinside_single, InternalBoundingBox, to_inside_index
import ...InternalBoundingBoxes: lower, upper
import ...Properties: GeometryLeafProperties
import ..Groups: GeometryTuple, GeometryVectorLike
import ..Repeats: Repeat
import ..Transformations: Transformation
import ..Transparents: Transparent

export SurfaceEstimate, estimate_volume, estimate_surface

contains_repeat(::Type{<:PhysicalGeometry}) = false
contains_repeat(::Type{<:Repeat}) = true
contains_repeat(::Type{<:GeometryVectorLike{N, P}}) where {N, P} = contains_repeat(P)
contains_repeat(::Type{<:GeometryTuple{N, P}}) where {N, P} =
    any(contains_repeat, P.parameters)
contains_repeat(::Type{<:Transformation{N, M, P}}) where {N, M, P} = contains_repeat(P)
contains_repeat(::Type{<:Transparent{N, P}}) where {N, P} = contains_repeat(P)

struct SurfaceEstimate{N, I <: Tuple}
    area::Float64
    positions::Vector{SVector{N, Float64}}
    normals::Vector{SVector{N, Float64}}
    indices::Vector{I}
    weights::Vector{Float64}
    density::Float64
    outer::Bool
end

function _inside_indices(geometry::PhysicalGeometry, position)
    has_single_inside(typeof(geometry)) ?
        (isinside_single(geometry, position) ? [()] : Tuple[]) :
        inside_indices(geometry, position)
end

function estimate_volume(
    geometry::PhysicalGeometry{N};
    bounding_box=nothing,
    nsamples::Integer=100_000,
    minimum_samples::Integer=1_000,
    maximum_attempts::Integer=5,
    rng::AbstractRNG=default_rng(),
) where {N}
    nsamples > 0 || throw(ArgumentError("nsamples must be positive"))
    minimum_samples >= 0 || throw(ArgumentError("minimum_samples must be nonnegative"))
    maximum_attempts > 0 || throw(ArgumentError("maximum_attempts must be positive"))
    has_inside(typeof(geometry)) || throw(ArgumentError("geometry has no inside volume"))
    contains_repeat(typeof(geometry)) &&
        throw(ArgumentError("Monte Carlo volume estimates do not support repeating geometries"))
    bounding_box = isnothing(bounding_box) ? InternalBoundingBox(geometry) : bounding_box

    box_size = upper(bounding_box) - lower(bounding_box)
    box_volume = prod(box_size)
    current_nsamples = nsamples
    inside_count = 0
    for attempt in 1:maximum_attempts
        inside_count = 0
        for _ in 1:current_nsamples
            position = SVector{N, Float64}(rand(rng, N)) .* box_size + lower(bounding_box)
            inside_count += !isempty(_inside_indices(geometry, position))
        end
        outside_count = current_nsamples - inside_count
        (inside_count >= minimum_samples && outside_count >= minimum_samples) && break
        current_nsamples *= 10
    end
    box_volume * inside_count / current_nsamples
end

function estimate_surface(
    geometry::PhysicalGeometry{N};
    bounding_box=nothing,
    density::Real=1.0,
    minimum_samples::Integer=1_000,
    maximum_attempts::Integer=5,
    outer::Bool=false,
) where {N}
    density > 0 || throw(ArgumentError("density must be positive"))
    minimum_samples >= 0 || throw(ArgumentError("minimum_samples must be nonnegative"))
    maximum_attempts > 0 || throw(ArgumentError("maximum_attempts must be positive"))
    contains_repeat(typeof(geometry)) &&
        throw(ArgumentError("Monte Carlo surface estimates do not support repeating geometries"))
    bounding_box = isnothing(bounding_box) ? InternalBoundingBox(geometry) : bounding_box

    current_density = Float64(density)
    positions = SVector{N, Float64}[]
    normals = SVector{N, Float64}[]
    indices = Tuple[]
    for _ in 1:maximum_attempts
        sampled_positions, sampled_indices = random_surface_positions(
            geometry,
            GeometryLeafProperties(current_density),
            bounding_box,
            1.0,
        )
        keep = trues(length(sampled_positions))
        if outer
            for index in eachindex(sampled_positions)
                full_index = sampled_indices[index]
                collision_indices = full_index[1:(end - 1)]
                current_inside = to_inside_index(geometry, collision_indices)
                keep[index] = all(
                    inside_index == current_inside
                    for inside_index in _inside_indices(geometry, sampled_positions[index])
                )
            end
        end

        accepted = findall(keep)
        positions = sampled_positions[accepted]
        indices = sampled_indices[accepted]
        normals = [
            get_intersection_params(
                geometry,
                position,
                position,
                (sampled_indices[index]..., 0.0),
            ).normal
            for (index, position) in zip(accepted, positions)
        ]
        length(positions) >= minimum_samples && break
        current_density *= 10
    end

    weights = fill(inv(current_density), length(positions))
    SurfaceEstimate(
        length(positions) * inv(current_density),
        positions,
        normals,
        indices,
        weights,
        current_density,
        outer,
    )
end

end
