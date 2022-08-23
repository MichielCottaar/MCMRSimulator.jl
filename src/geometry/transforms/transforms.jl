"""
Transformations of the [`BaseObstruction`](@ref).

In order the following transformations are applied to each obstruction:
- [`Shifted{N}`](@ref): Applies shift in relevant dimensions for [`BaseObstruction`](@ref).
- [`Repeated{N}`](@ref): Repeat obstruction infinitely many times.
- [`Transformed{N}`](@ref): applies rotation & shift before projecting to lower dimensionality.
"""
abstract type TransformObstruction{N} <: Obstruction{N} end


include("repeat.jl")
include("shift.jl")
include("transform.jl")