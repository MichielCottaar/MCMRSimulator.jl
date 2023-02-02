"""
Supertype of any obstruction to the free diffusion of water. Some might also generate off-resonance fields.

There are two types of obstructions:
- [`BaseObstruction`](@ref) with the basic obstructions (walls, sphere, cylinders, meshes)
- [`TransformObstruction`](@ref), which transform the base obstructions
"""
abstract type Obstruction{N} end

include("movement.jl")  # Bounding boxes
include("bounding_box.jl")  # Bounding boxes
include("properties.jl")    # Obstruction properties
include("collision.jl")     # Collision object returned by detect_collision
include("grid.jl")          # Repeating grids
include("base/base.jl")     # Base obstructions
include("transform.jl")     # Transformations of base obstructions
include("geometry_struct.jl")     # Transformations of base obstructions
include("diffuse.jl")       # Actual diffusion algorithm


"""
    project(position, transform::TransformObstruction)

Computes the position in the space of the obstructions wrapped by the [`TransformObstruction`](@ref).

    project(position, grid::GridShape)

Computes the voxel index for the position on the [`GridShape`](@ref). 
This will return a result even if the point is outside of the grid. Use [`isinside`](@ref)(position, grid) to check that.
"""
function project end