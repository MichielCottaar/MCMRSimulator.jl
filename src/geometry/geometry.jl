"""
Supertype of any obstruction to the free diffusion of water. Some might also generate off-resonance fields.

There are two types of obstructions:
- [`BaseObstruction`](@ref) with the basic obstructions (walls, sphere, cylinders, meshes)
- [`TransformObstruction`](@ref), which transform the base obstructions
"""
abstract type Obstruction{N} end

include("movement.jl")  # Bounding boxes
include("bounding_box.jl")  # Bounding boxes
include("collision.jl")     # Collision object returned by detect_collision
include("reflection.jl")
include("grid.jl")          # Repeating grids
include("base/base.jl")     # Base obstructions
include("transform.jl")     # Transformations of base obstructions
include("geometry_struct.jl")     # Transformations of base obstructions


