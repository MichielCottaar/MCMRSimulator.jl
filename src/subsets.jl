"""
Support for selecting a subset of the total signal.

Types:
- [`Subset`](@ref)

Functions:
- [`get_subset`](@ref)
"""
module Subsets
import Base: @kwdef
import ..Spins: Snapshot, Spin, stuck_to
import ..Reflections: previous_hit
import ..Simulations: Simulation
import ..Geometries: fix
import ..Geometries.Internal: FixedGeometry, distance_to_surface, isinside

_arguments = """
- `geometry`: user geometry to use for inside and bound tests (default: the simulation geometry).
- `bound`: set to true to return only bound spins, to false to return only free spins (default: whether spins are bound is not relevant).
- `inside`: set to true to return only spins inside the geometry, to false to return only spins outside of the geometry (default: whether spins are inside or outside is not relevant). It can also be set to a positive integer number. Only spins that are inside that exact number of obstructions will be returned.
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

get_subset(snapshot::Snapshot, simulation::Union{Simulation, FixedGeometry}; kwargs...) =
    get_subset(snapshot, simulation, Subset(; kwargs...))

function get_subset(snapshot::Snapshot, simulation::Simulation, subset::Subset)
    if isnothing(subset.geometry)
        return _get_subset(snapshot, simulation.geometry, subset, false)
    end
    _get_subset(snapshot, fix(subset.geometry), subset, true)
end

get_subset(snapshot::Snapshot, geometry::FixedGeometry, subset::Subset) =
    _get_subset(snapshot, geometry, subset, false)

function _get_subset(snapshot::Snapshot{N}, geometry::FixedGeometry, subset::Subset, explicit_geometry::Bool) where {N}
    if isnothing(subset.bound) && isnothing(subset.inside)
        return snapshot
    end

    include = ones(Bool, length(snapshot))
    if !isnothing(subset.bound)
        function _isbound(spin::Spin)
            if explicit_geometry
                return !isempty(stuck_to(spin)) && distance_to_surface(geometry, spin.position) <= sqrt(eps())
            end
            return !isempty(stuck_to(spin))
        end
        include .&= (_isbound.(snapshot) .== subset.bound)
    end
    if !isnothing(subset.inside)
        function _number_isinside(spin::Spin)
            length(isinside(geometry, spin.position, previous_hit(spin.reflection)))
        end
        if subset.inside isa Bool
            include .&= ((_number_isinside.(snapshot) .> 0) .== subset.inside)
        else
            include .&= (_number_isinside.(snapshot) .== subset.inside)
        end
    end
    return snapshot[include]
end


end
