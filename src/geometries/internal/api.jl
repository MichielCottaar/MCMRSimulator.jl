"""Public operations supported by the internal geometry engine."""

import StaticArrays: SVector
import ..BoundingBoxes: BoundingBox
import .Indices: ObstructionIndex
import .InternalBoundingBoxes: InternalBoundingBox
import .InternalBoundingBoxes
import .PhysicalGeometries: PhysicalGeometry, Intersection, flip, detect_intersection, random_surface_positions, inside_indices, size_scale, geometry_mesh
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
    IsInside

Identifies which obstructions this position/spin is inside of.
Use `Base.length` to get the number of obstructions.
"""
struct IsInside
    inside_of::Vector{ObstructionIndex}
end

Base.length(ii::IsInside) = length(ii.inside_of)


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
isinside(geometry::FixedGeometry, position::SVector{3, Float64}, intersection::Intersection=Intersection{3}()) =
    IsInside(inside_indices_for_any_type(geometry.geometry, position, intersection))

"""
    detect_intersection(geometry, start, destination, previous_intersection)

Return the first intersection between the path from `start` to `destination`
and `geometry`.
"""
function detect_intersection(
    geometry::FixedGeometry,
    start::SVector{3, Float64},
    destination::SVector{3, Float64},
    previous_intersection::Intersection=Intersection{3}(),
)
    detect_intersection(geometry.geometry, start, destination, previous_intersection)
end

function _filter_to_box(positions, values, bounding_box)
    keep = [all(position .>= InternalBoundingBoxes.lower(bounding_box)) &&
        all(position .<= InternalBoundingBoxes.upper(bounding_box)) for position in positions]
    positions[keep], values[keep]
end

function random_surface_positions(
    geometry::FixedGeometry,
    bounding_box::InternalBoundingBox{3},
    volume_density::Number,
)
    isnothing(geometry.surface) && return SVector{3, Float64}[], Intersection{3}[]
    positions, intersections = random_surface_positions(
        geometry.geometry, geometry.surface.density, bounding_box, volume_density)
    _filter_to_box(positions, intersections, bounding_box)
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
    isnothing(geometry.surface) && return Inf

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
    isnothing(geometry.surface) && return 0.
    values = (value for value in all_property_values(geometry.surface.permeability) if !isinf(value))
    maximum(values; init=0.)
end

"""Return the largest surface-relaxation value in `geometry`."""
function max_surface_relaxation(geometry::FixedGeometry)
    isnothing(geometry.surface) && return 0.
    maximum(all_property_values(geometry.surface.surface_relaxation); init=0.)
end

"""Return the smallest dwell time in `geometry` when a bound pool exists."""
function min_dwell_time(geometry::FixedGeometry)
    isnothing(geometry.surface) && return Inf

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
    @assert !Base.isempty(intersection) "permeability requires a non-empty intersection"
    intersection.hit_gap && return Inf
    get_value(geometry.surface.permeability, intersection.obstruction_index)
end

"""Return the surface-relaxation value associated with a collision state."""
function surface_relaxation(geometry::FixedGeometry, intersection::Intersection)
    @assert !Base.isempty(intersection) "surface_relaxation requires a non-empty intersection"
    intersection.hit_gap && return 0.
    get_value(geometry.surface.surface_relaxation, intersection.obstruction_index)
end

"""Return the surface-density value associated with a collision state."""
function surface_density(geometry::FixedGeometry, intersection::Intersection)
    @assert !Base.isempty(intersection) "surface_density requires a non-empty intersection"
    intersection.hit_gap && return 0.
    get_value(geometry.surface.density, intersection.obstruction_index)
end

"""Return the dwell-time value associated with a collision state."""
function dwell_time(geometry::FixedGeometry, intersection::Intersection)
    @assert !Base.isempty(intersection) "dwell_time requires a non-empty intersection"
    intersection.hit_gap && return 0.
    get_value(geometry.surface.dwell_time, intersection.obstruction_index)
end

# Generate the shared volume-plus-surface lookup for MRI properties.
for symbol in (:R1, :R2, :off_resonance)
    @eval begin
        function $symbol(
            geometry::FixedGeometry,
            inside::IsInside,
            intersection::Intersection=Intersection{3}(),
        )
            volume = get_value(geometry.volume.$symbol, inside.inside_of)
            surface = Base.isempty(intersection) ?
                0.0 :
                get_value(geometry.surface.$symbol, intersection.obstruction_index)
            volume + surface
        end

        function $symbol(
            geometry::FixedGeometry,
            position::SVector{3, Float64},
            intersection::Intersection=Intersection{3}(),
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
