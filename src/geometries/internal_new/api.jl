"""Methods and types implemented by the internal geometry engine."""
module API

export FixedGeometry, CollisionState, BoundingBox,
    fix, isinside, detect_intersection, random_surface_positions, geometry_mesh,
    size_scale, max_timestep_sticking, max_permeability_non_inf,
    max_surface_relaxation, min_dwell_time,
    permeability, surface_relaxation, surface_density, dwell_time,
    mri_properties, R1, R2, off_resonance,
    susceptibility_off_resonance, off_resonance_gradient,
    has_intersection, previous_hit,
    direction, reflect

"""
    FixedGeometry

Opaque fixed geometry state used by the simulator. Implementations may store
collision, MRI-property, susceptibility, and plotting data in any way they
choose. Callers should use only the functions defined in this module.
"""
abstract type FixedGeometry end

"""
    CollisionState

Opaque result of an intersection query or state describing a reflected or
surface-bound spin. Use `Base.isempty` to determine whether the state contains
an intersection; an empty state represents no collision.
"""
abstract type CollisionState end

"""
    fix(geometry; kwargs...)

Convert user-facing geometry into one `FixedGeometry`. The result must include
all state needed by collision detection, MRI-property queries, susceptibility
queries, timestep limits, and plotting.
"""
function fix end

"""
    isinside(geometry, position, reflection)

Return the obstructions containing `position`. The result must identify each
hit sufficiently for the simulator to apply volume properties. An omitted
reflection is interpreted as an unbound spin by the implementation.
"""
function isinside end

"""
    detect_intersection(geometry, start, destination, previous_hit)

Find the first intersection of the path from `start` to `destination`.
Return a `CollisionState`. Use `Base.isempty` on the result to determine
whether there was no hit.
"""
function detect_intersection end

"""
    random_surface_positions(geometry, bounding_box, volume_density)

Sample surface positions in `bounding_box` at the requested volume density.
Each result contains a position, inward unit normal, geometry index, and
obstruction index.
"""
function random_surface_positions end

"""
    geometry_mesh(geometry)

Return render-ready mesh data containing vertices and triangles. Plotting code
must not inspect the internal geometry representation.
"""
function geometry_mesh end

"""Return the smallest relevant obstruction size scale in `geometry`."""
function size_scale end

"""Return the maximum timestep allowed by surface sticking constraints."""
function max_timestep_sticking end

"""Return the largest finite permeability in `geometry`, or zero if none exists."""
function max_permeability_non_inf end

"""Return the largest surface-relaxation value in `geometry`."""
function max_surface_relaxation end

"""Return the smallest dwell time for a geometry with a bound pool."""
function min_dwell_time end

"""Return the permeability at a `CollisionState`."""
function permeability end

"""Return the surface-relaxation value at a `CollisionState`."""
function surface_relaxation end

"""Return the surface-density value at a `CollisionState`."""
function surface_density end

"""Return the dwell-time value at a `CollisionState`."""
function dwell_time end

"""
    mri_properties(geometry, global_properties, position, reflection)

Return an object with `R1`, `R2`, and `off_resonance` properties for a position
and optional surface reflection. Callers should not depend on its type.
"""
function mri_properties end

"""Return longitudinal relaxation at `position` in `geometry`."""
function R1 end

"""Return transverse relaxation at `position` in `geometry`."""
function R2 end

"""Return non-susceptibility off-resonance at `position` in `geometry`."""
function off_resonance end

"""Return susceptibility-induced off-resonance at a position or along a path."""
function susceptibility_off_resonance end

"""Return the maximum susceptibility off-resonance gradient for `geometry`."""
function off_resonance_gradient end

"""Return whether a `CollisionState` represents a real collision."""
function has_intersection end

"""Return the previous-hit information needed by `detect_intersection`."""
function previous_hit end

"""Return the remaining displacement for a `CollisionState` over a time interval."""
function direction end

"""
    reflect(collision, direction, ratio_displaced, time_moved, distance_moved;
            permeable=false)

Create a `CollisionState` after a collision. A permeable collision continues
through the surface; otherwise the direction is reflected.
"""
function reflect end

end
