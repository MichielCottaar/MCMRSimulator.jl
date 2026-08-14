"""Public operations supported by the internal geometry engine."""

import .Indices: ObstructionIndex
import .PhysicalGeometries: PhysicalGeometry, Intersection, surface_sampling
import .Properties: all_property_values
import ...Properties: stick_probability

export FixedGeometry, Intersection, ObstructionIndex,
    isinside, detect_intersection, random_surface_positions, surface_sampling, geometry_mesh,
    size_scale, max_timestep_sticking, max_permeability_non_inf,
    max_surface_relaxation, min_dwell_time,
    permeability, surface_relaxation, surface_density, dwell_time,
    mri_properties, R1, R2, off_resonance,
    susceptibility_off_resonance, off_resonance_gradient,
    direction, reflect

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

"""
    isinside(geometry, position[, reflection])

Return the obstruction indices containing `position`. A reflection may be
provided when the particle is already attached to a surface.
"""
function isinside end

"""
    detect_intersection(geometry, start, destination, previous_intersection)

Return the first intersection between the path from `start` to `destination`
and `geometry`.
"""
function detect_intersection end

"""
    random_surface_positions(geometry, bounding_box, volume_density)

Sample surface positions in `bounding_box` at the requested volume density.
"""
function random_surface_positions end

"""
    geometry_mesh(geometry)

Return render-ready mesh data for `geometry`.
"""
function geometry_mesh end

"""Return the smallest relevant obstruction size in `geometry`."""
function size_scale end

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
function permeability end

"""Return the surface-relaxation value associated with a collision state."""
function surface_relaxation end

"""Return the surface-density value associated with a collision state."""
function surface_density end

"""Return the dwell-time value associated with a collision state."""
function dwell_time end

"""
    mri_properties(geometry, global_properties, position, reflection)

Return the `R1`, `R2`, and off-resonance values at `position`, optionally taking
a surface reflection into account.
"""
function mri_properties end

"""Return longitudinal relaxation at `position` in `geometry`."""
function R1 end

"""Return transverse relaxation at `position` in `geometry`."""
function R2 end

"""Return non-susceptibility off-resonance at `position` in `geometry`."""
function off_resonance end

"""Return susceptibility-induced off-resonance for `geometry`."""
function susceptibility_off_resonance end

"""Return the maximum susceptibility off-resonance gradient in `geometry`."""
function off_resonance_gradient end

"""Return the remaining displacement associated with a collision state."""
function direction end

"""
    reflect(collision, direction, ratio_displaced, time_moved, distance_moved;
            permeable=false)

Return the collision state after reflecting or passing through a surface.
"""
function reflect end
