"""
Methods to extract properties from the [`FixedGeometry`](@ref).

MRI properties (using [`MRIProperties`](@ref)):
- [`R1`](@ref)
- [`R2`](@ref)
- [`off_resonance`](@ref)

Collision properties:
- [`permeability`](@ref)
- [`surface_relaxation`](@ref)
- [`surface_density`](@ref)
- [`dwell_time`](@ref)

Others:
- [`max_timestep_sticking`](@ref)
"""
module Properties

import StaticArrays: SVector
import ....Properties: GlobalProperties, R1, R2, off_resonance, stick_probability
import ..FixedObstructionGroups: FixedGeometry, FixedObstructionGroup, isinside
import ..Intersections: Intersection
import ..Reflections: Reflection, empty_reflection

"""
    MRIProperties(R1, R2, off_resonance)
    MRIProperties(full_geometry, inside_geometry, global_properties, position, reflection)

Determines the R1, R2, and off-resonance of a particle at given `position` and bound at [`Reflection`](@ref).
"""
mutable struct MRIProperties
    R1 :: Float64
    R2 :: Float64
    off_resonance :: Float64
end

function MRIProperties(full_geometry::FixedGeometry, inside_geometry::FixedGeometry, glob::GlobalProperties, position::SVector{3, Float64}, stuck_to::Reflection)
    res = MRIProperties(glob.R1, glob.R2, glob.off_resonance)

    # first check stuck properties
    if ~iszero(stuck_to.geometry_index)
        group = full_geometry[stuck_to.geometry_index]
        res.R1 += get_value(Val(:R1), group.surface, stuck_to.obstruction_index)
        res.R2 += get_value(Val(:R2), group.surface, stuck_to.obstruction_index)
        res.off_resonance += get_value(Val(:off_resonance), group.surface, stuck_to.obstruction_index)
    end

    for group in inside_geometry
        props = group.volume
        indices = isinside(group, position, stuck_to)
        for index in indices
            res.R1 += get_value(Val(:R1), props, index)
            res.R2 += get_value(Val(:R2), props, index)
            res.off_resonance += get_value(Val(:off_resonance), props, index)
        end
    end
    return res
end

function get_value(::Val{S}, properties::NamedTuple, index::Integer) where {S}
    res = getproperty(properties, S)
    if res isa Vector
        return res[index]
    else
        return res
    end
end

for symbol in (:R1, :R2, :off_resonance)
    @eval begin
        function $symbol(position::AbstractVector, geometry::FixedGeometry, properties::GlobalProperties=GlobalProperties(), stuck_to::Reflection=empty_reflection)
            MRIProperties(geometry, geometry, properties, SVector{3, Float64}(position), stuck_to).$symbol
        end
    end
end


for symbol in (:permeability, :surface_relaxation, :dwell_time, :surface_density)
    @eval begin
        """
            $($symbol)(geometry, reflection)
        
        Returns the $($symbol) experienced by the spin hitting the surface represented by a [`Reflection`](@ref).
        The `geometry` has to be a [`FixedGeometry`](@ref).
        """
        function $symbol(group::FixedObstructionGroup, has_hit::Int)
            value = group.surface.$symbol
            if value isa Vector
                return value[has_hit]
            else
                return value
            end
        end

        function $symbol(geometry::FixedGeometry, has_hit::Tuple{Int, Int})
            $symbol(geometry[has_hit[1]], has_hit[2])
        end

        function $symbol(geometry::FixedGeometry, has_hit::Union{Intersection, Reflection})
            $symbol(geometry, (has_hit.geometry_index, has_hit.obstruction_index))
        end
    end
end

"""
    max_timestep_sticking(geometry, diffusivity)

Returns the maximum timestep that can be used while:
1. keeping [`stick_probability`](@ref) lower than 25%
2. and keeping the timestep below the [`dwell_time`](@ref)
"""
function max_timestep_sticking(geometry::FixedGeometry, diffusivity::Number)
    minimum([max_timestep_sticking(group, diffusivity) for group in geometry])
end

max_timestep_sticking(geometry::FixedGeometry{0}, diffusivity::Number) = Inf

function max_timestep_sticking(group::FixedObstructionGroup, diffusivity::Number)
    get_dwell_time(surface_density, dwell_time) = iszero(surface_density) ? Inf : dwell_time
    return min(
        maximum(stick_probability.(group.surface.surface_density, group.surface.dwell_time, diffusivity, 1))^-2 / 4,
        minimum(get_dwell_time.(group.surface.surface_density, group.surface.dwell_time))
    )
end

"""
    max_permeability_non_inf(geometry)

Returns the largest permeability within the geometry that is not infinite.
"""
function max_permeability_non_inf(geometry::FixedGeometry)
    maximum(max_permeability_non_inf.(geometry))
end

max_permeability_non_inf(geometry::FixedGeometry{0}) = 0.

function max_permeability_non_inf(group::FixedObstructionGroup)
    if group.surface.permeability isa AbstractVector
        return filter(a -> ~isinf(a), group.surface.permeability)
    else
        return isinf(group.surface.permeability) ? 0. : group.surface.permeability
    end
end

end