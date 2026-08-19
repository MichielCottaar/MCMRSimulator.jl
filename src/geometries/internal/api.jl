"""Public operations supported by the internal geometry engine."""

import StaticArrays: SVector
import ..BoundingBoxes: BoundingBox
import .InternalBoundingBoxes: InternalBoundingBox
import .InternalBoundingBoxes
import .PhysicalGeometries: PhysicalGeometry, find_intersection, get_intersection_params, random_surface_positions, inside_indices, size_scale, geometry_mesh, to_property_index, inside_indices_eltype
import .PhysicalGeometries.Groups: GeometryTuple, inside_indices_for_any_type
import .PhysicalGeometries.Transparents: SizeScaleOverride
import .Properties: all_property_values, get_value
import .Susceptibility: susceptibility_off_resonance, off_resonance_gradient
import ...Properties: stick_probability
import .RayGridIntersection: ray_grid_intersections

export FixedGeometry, Intersection, flip, IsInside, collision_normal,
    isinside, detect_intersection, random_surface_positions, geometry_mesh,
    ray_grid_intersections,
    size_scale, SizeScaleOverride, max_timestep_sticking, max_permeability_non_inf,
    max_surface_relaxation, min_dwell_time,
    permeability, surface_relaxation, surface_density, dwell_time,
    R1, R2, off_resonance,
    susceptibility_off_resonance, off_resonance_gradient

"""
    FixedGeometry

The fixed representation of a user-defined geometry used by the simulation.
It contains the geometry and all fixed state needed for geometry queries.
"""
struct FixedGeometry{G <: PhysicalGeometry{3}, V, S, O}
    geometry::G
    volume::V
    surface::S
    susceptibility::O
end

# A fixed geometry represents one user geometry unless it is a tuple of groups.
Base.length(geometry::FixedGeometry) =
    geometry.geometry isa GeometryTuple ? length(geometry.geometry) : 1


"""
    Intersection(distance, indices, property_indices, normal, inside, hit_gap)

Intersection between a path and a physical geometry. 

- `distance` is normalized to the path from start (0) to destination (1) 
- `indices` is a tuple containing any information required to find the obstruction of the intersection again, excluding `inside` and `distance`
- `normal` is the normal of the surface hit by the collision
- `inside` is whether the surface got hit on the "inside" of the obstruction (arbitrarily defined for some obstructions such as infinite walls)
- `hit_gap` is whether the collision actually hit something that represents a gap in the obstruction. In that case the spin will be left through unaltered except for an update to whether it is inside/outside of that obstruction.
"""
struct Intersection{T, PT}
    distance::Float64
    indices::T
    property_indices::PT
    normal::SVector{3, Float64}
    inside::Bool
    hit_gap::Bool
end

"""
    flip(intersection)

Flip an intersection to the opposite side of its surface while preserving its
distance, obstruction indices, property indices, and gap state.
"""
function flip(intersection::Intersection{T, PT}) where {T, PT}
    Intersection{T, PT}(
        intersection.distance,
        intersection.indices,
        intersection.property_indices,
        -intersection.normal,
        !intersection.inside,
        intersection.hit_gap,
    )
end


"""
    detect_intersection(fixed_geometry, start, dest, previous_intersection=nothing) -> Optional[Intersection]

Finds the closest intersection between `start` and `dest` with `fixed_geometry`.

If no intersection is found, `nothing` is returned instead.
`previous_intersection` represents the intersection that just finished, which should not be immediately returned again.
"""
function detect_intersection(fixed_geometry::FixedGeometry, start::SVector{3, Float64}, dest::SVector{3, Float64}, previous_intersection=nothing)
    full_indices = find_intersection(
        fixed_geometry.geometry,
        start,
        dest,
        isnothing(previous_intersection) ?
        nothing :
        (previous_intersection.indices..., previous_intersection.inside, previous_intersection.distance),
    )
    if isnothing(full_indices)
        return nothing
    end
    params = get_intersection_params(fixed_geometry.geometry, start, dest, full_indices)
    inside = full_indices[end - 1]
    distance = full_indices[end]
    indices = full_indices[1:end - 2]
    property_indices = to_property_index(fixed_geometry.geometry, indices)
    return Intersection{typeof(indices), typeof(property_indices)}(
        distance,
        indices,
        property_indices,
        params.normal,
        inside,
        params.hit_gap
    )
end



"""
    IsInside

Identifies which obstructions this position/spin is inside of.
Use `Base.length` to get the number of obstructions.
"""
struct IsInside{N, T}
    inside_of::SVector{N, T}
end

Base.length(::IsInside{N}) where {N} = N


"""
    collision_normal(intersection)

Returns the normal of the intersection.
"""
collision_normal(intersection::Intersection) = intersection.normal

"""
    isinside(geometry, position[, intersection])

Return the obstruction indices containing `position`. 

An intersection may be provided when the particle is already attached or is currently hitting a surface.
"""
function isinside(
    geometry::FixedGeometry,
    position::SVector{3, Float64},
    intersection=nothing,
)
    as_vec = inside_indices_for_any_type(geometry.geometry, position, intersection)
    expected_type = inside_indices_eltype(typeof(geometry.geometry))
    IsInside(SVector{length(as_vec), expected_type}(as_vec))
end

isinside(geometry::FixedGeometry, position::SVector{3, Float64}, intersection::Intersection) =
    isinside(geometry, position, (intersection.indices..., intersection.inside, intersection.distance))


function _filter_to_box(positions, values, bounding_box)
    keep = [all(position .>= InternalBoundingBoxes.lower(bounding_box)) &&
        all(position .<= InternalBoundingBoxes.upper(bounding_box)) for position in positions]
    positions[keep], values[keep]
end

function random_surface_positions(
    fixed_geometry::FixedGeometry,
    bounding_box::InternalBoundingBox{3},
    volume_density::Number,
)
    all_positions, all_indices = random_surface_positions(
        fixed_geometry.geometry, fixed_geometry.surface.density, bounding_box, volume_density
    )
    
    positions, indices = _filter_to_box(all_positions, all_indices, bounding_box)
    intersections = map(positions, indices) do pos, full_index
        inside = full_index[end]
        index = full_index[1:end - 1]
        params = get_intersection_params(fixed_geometry.geometry, pos, pos, (full_index..., 0.))
        property_index = to_property_index(fixed_geometry.geometry, index)
        return Intersection{typeof(index), typeof(property_index)}(
            0.,
            index,
            property_index,
            params.normal,
            inside,
            params.hit_gap
        )
    end
    return (positions, intersections)
end

function random_surface_positions(
    geometry::FixedGeometry,
    bounding_box::BoundingBox,
    volume_density::Number,
)
    random_surface_positions(geometry, InternalBoundingBox(bounding_box), volume_density)
end

"""
    geometry_mesh(geometry)

Return render-ready mesh data for `geometry`.
"""
function geometry_mesh(geometry::FixedGeometry; bounding_box=nothing, kwargs...)
    bounding_box = bounding_box isa BoundingBox ? InternalBoundingBox(bounding_box) : bounding_box
    geometry_mesh(geometry.geometry; bounding_box, kwargs...)
end

function geometry_mesh(geometry::FixedGeometry, bounding_box::BoundingBox; kwargs...)
    geometry_mesh(geometry; bounding_box=InternalBoundingBox(bounding_box), kwargs...)
end

function geometry_mesh(geometry::PhysicalGeometry, bounding_box::BoundingBox; kwargs...)
    geometry_mesh(geometry; bounding_box=InternalBoundingBox(bounding_box), kwargs...)
end

"""Return the smallest relevant obstruction size in `geometry`."""
size_scale(geometry::FixedGeometry) = size_scale(geometry.geometry)

"""Return the largest timestep permitted by surface-sticking constraints."""
function max_timestep_sticking(geometry::FixedGeometry, diffusivity::Number, scaling)
    log_probabilities = [
        log(1 - stick_probability(density, dwell, diffusivity, 1))
        for (density, dwell) in all_property_values(
            geometry.surface.density,
            geometry.surface.dwell_time,
        )
        if !iszero(density)
    ]
    isempty(log_probabilities) && return Inf
    scaling * maximum(log_probabilities)^(-2)
end

"""Return the largest finite permeability in `geometry`, or zero if absent."""
function max_permeability_non_inf(geometry::FixedGeometry)
    values = (value for value in all_property_values(geometry.surface.permeability) if !isinf(value))
    maximum(values; init=0.)
end

"""Return the largest surface-relaxation value in `geometry`."""
function max_surface_relaxation(geometry::FixedGeometry)
    maximum(all_property_values(geometry.surface.surface_relaxation); init=0.)
end

"""Return the smallest dwell time in `geometry` when a bound pool exists."""
function min_dwell_time(geometry::FixedGeometry)
    dwell_times = (
        dwell
        for (density, dwell) in all_property_values(
            geometry.surface.density,
            geometry.surface.dwell_time,
        )
        if !iszero(density)
    )
    minimum(dwell_times; init=Inf)
end

"""Return the permeability associated with a collision state."""
function permeability(geometry::FixedGeometry, intersection::Intersection)
    intersection.hit_gap && return Inf
    get_value(geometry.surface.permeability, intersection.property_indices)
end

"""Return the surface-relaxation value associated with a collision state."""
function surface_relaxation(geometry::FixedGeometry, intersection::Intersection)
    intersection.hit_gap && return 0.
    get_value(geometry.surface.surface_relaxation, intersection.property_indices)
end

"""Return the surface-density value associated with a collision state."""
function surface_density(geometry::FixedGeometry, intersection::Intersection)
    intersection.hit_gap && return 0.
    get_value(geometry.surface.density, intersection.property_indices)
end

"""Return the dwell-time value associated with a collision state."""
function dwell_time(geometry::FixedGeometry, intersection::Intersection)
    intersection.hit_gap && return 0.
    get_value(geometry.surface.dwell_time, intersection.property_indices)
end

# Generate the shared volume-plus-surface lookup for MRI properties.
for symbol in (:R1, :R2, :off_resonance)
    @eval begin
        function $symbol(
            geometry::FixedGeometry,
            inside::IsInside,
            intersection=nothing,
        )
            volume = get_value(geometry.volume.$symbol, inside.inside_of)
            surface = isnothing(intersection) ?
                0.0 :
                get_value(geometry.surface.$symbol, intersection.property_indices)
            volume + surface
        end

        function $symbol(
            geometry::FixedGeometry,
            position::SVector{3, Float64},
            intersection=nothing,
        )
            $symbol(geometry, isinside(geometry, position, intersection), intersection)
        end
    end
end

"""Return susceptibility-induced off-resonance for `geometry`."""
function susceptibility_off_resonance(
    geometry::FixedGeometry,
    position::SVector{3, Float64},
    inside::Union{Nothing, Bool}=nothing,
)
    susceptibility_off_resonance(geometry.susceptibility, position, inside)
end

"""Return the maximum susceptibility off-resonance gradient in `geometry`."""
function off_resonance_gradient(geometry::FixedGeometry, B0)
    off_resonance_gradient(geometry.susceptibility, B0)
end
