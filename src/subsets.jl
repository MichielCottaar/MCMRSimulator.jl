"""
Support for selecting a subset of the total signal.

Types:
- [`Subset`](@ref)

Functions:
- [`get_subset`](@ref)
"""
module Subsets
import Base: @kwdef
import ..Spins: Snapshot, Spin
import ..Simulations: Simulation, isinside, stuck

_arguments = """
- `geometry`: user geometry or view selecting obstruction IDs in the simulation geometry (default: the full simulation geometry). IDs that do not match the simulation geometry do not contribute to the result.
- `bound`: set to true to return only spins bound to a selected obstruction, to false to return spins not bound to one (default: whether spins are bound is not relevant).
- `inside`: set to true to return spins with at least one matching inside index, or false to return spins with none (default: whether spins are inside is not relevant). An integer selects spins with exactly that many matching inside indices.
"""

"""
    Subset(; geometry=nothing, bound=nothing, inside=nothing)

This creates a helper object to extract a subset of a [`Snapshot`](@ref) from the total snapshot.
It defines which spins should be kept.
This definition is determined by:
$(_arguments)
"""
@kwdef struct Subset
    geometry = nothing
    bound :: Union{Nothing, Bool} = nothing
    inside :: Union{Nothing, Bool, Int} = nothing
end


"""
    get_subset(snapshot, simulation, subset)
    get_subset(snapshot, simulation; geometry=nothing, bound=nothing, inside=nothing)

Returns a subset of the [`Snapshot`](@ref) from the [`Simulation`](@ref) that obeys some specific properties.
These properties can be either defined by a [`Subset`](@ref) object or a set of keyword arguments.

These keyword arguments are:
$(_arguments)
"""
get_subset(spins::AbstractVector{<:Spin}, args...; kwargs...) = get_subset(Snapshot(spins), args...; kwargs...).spins

get_subset(snapshot::Snapshot, simulation::Simulation; kwargs...) =
    get_subset(snapshot, simulation, Subset(; kwargs...))

function get_subset(snapshot::Snapshot, simulation::Simulation, subset::Subset)
    if isnothing(subset.bound) && isnothing(subset.inside)
        return snapshot
    end

    include = ones(Bool, length(snapshot))
    if !isnothing(subset.bound)
        bound = isnothing(subset.geometry) ?
            stuck(simulation, snapshot) :
            stuck(simulation, subset.geometry, snapshot)
        include .&= (bound .== subset.bound)
    end
    if !isnothing(subset.inside)
        number_inside = isnothing(subset.geometry) ?
            isinside(simulation, snapshot) :
            isinside(simulation, subset.geometry, snapshot)
        if subset.inside isa Bool
            include .&= ((number_inside .> 0) .== subset.inside)
        else
            include .&= (number_inside .== subset.inside)
        end
    end
    return snapshot[include]
end


end
