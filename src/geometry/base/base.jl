"""Base obstruction type hindering free diffusion.

These obstructions are always perfectly aligned with the main axes and centered at the origin (except for [`Mesh`](@ref)).
They can be moved/rotated/repeated by applying [`TransformObstruction`](@ref).

The dimensionality N indicates the dimensionality of the input data 
(1 for [`Wall`](@ref), 2 for [`Cylinder`](@ref), 3 for [`Sphere`](@ref) or [`Mesh`](@ref)).
"""
abstract type BaseObstruction{N} <: Obstruction{N} end

ObstructionProperties(obstruction :: BaseObstruction) = obstruction.properties

Base.copy(o::BaseObstruction) = @set o.properties.id = uuid1()

include("random_draws.jl")
include("cylinder.jl")
include("sphere.jl")
include("wall.jl")
include("mesh.jl")
include("annulus.jl")