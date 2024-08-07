"""
Internal representation of the tissue geometry.

All objects are immutable and optimised for both memory usage and computational speed.
"""
module Internal
include("ray_grid_intersection.jl")
include("bounding_boxes.jl")
include("obstructions/obstructions.jl")
include("hit_grids.jl")
include("intersections.jl")
include("reflections.jl")
include("fixed_obstruction_groups.jl")
include("properties.jl")
include("susceptibility/susceptibility.jl")
include("isinside_mesh.jl")

import .FixedObstructionGroups: FixedObstructionGroup, FixedGeometry, FixedMesh, repeating, isinside, detect_intersection
import .Reflections: Reflection, direction, previous_hit, has_hit, has_intersection, empty_reflection
import .Intersections: Intersection, empty_intersection
import .HitGrids: HitGrid, detect_intersection_grid
import .Obstructions: 
    FixedObstruction, ObstructionIntersection, empty_obstruction_intersections,
    Wall, Round, Cylinder, Sphere, Triangle, IndexTriangle, FullTriangle, Shift,
    has_inside, radius, obstruction_type, size_scale, get_quadrant, random_surface_positions
import .BoundingBoxes: BoundingBox, could_intersect, lower, upper
import .RayGridIntersection: ray_grid_intersections
import .Properties: R1, R2, off_resonance, permeability, surface_relaxivity, surface_density, dwell_time, max_timestep_sticking, MRIProperties
import .Susceptibility: FixedSusceptibility, susceptibility_off_resonance, off_resonance_gradient
import .IsInsideMesh: isinside_grid

end