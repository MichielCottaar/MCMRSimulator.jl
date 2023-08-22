"""
Methods to extract properties from the [`FixedGeometry`](@ref).

Using [`MRIProperties`](@ref):
- [`R1`](@ref)
- [`R2`](@ref)
- [`off_resonance`](@ref)

Using [`collision_property`](@ref):
- [`permeability`](@ref)
- [`surface_relaxivity`](@ref)
- [`surface_density`](@ref)
- [`dwell_time`](@ref)

Others:
- [`max_timestep_sticking`](@ref)
"""
module Properties

import StaticArrays: SVector
import ....Properties: GlobalProperties, R1, R2, off_resonance, permeability, surface_relaxivity, surface_density, dwell_time, stick_probability
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
    res = MRIProperties(NaN, NaN, NaN)

    finished() = ~(isnan(res.R1) || isnan(res.R2) || isnan(res.off_resonance))

    # first check stuck properties
    if ~iszero(stuck_to.geometry_index)
        update_stuck_properties!(res, full_geometry[stuck_to.geometry_index], stuck_to.obstruction_index)
        if finished()
            return res
        end
    end

    for group in inside_geometry
        if (
            (isnan(res.R1) && hasproperty(group.volume, :R1)) ||
            (isnan(res.R2) && hasproperty(group.volume, :R2)) ||
            (isnan(res.off_resonance) && hasproperty(group.volume, :off_resonance))
        )
            update_inside_properties!(res, group, position, group.parent_index == stuck_to.geometry_index ? stuck_to.obstruction_index : empty_reflection)
        end
        if finished()
            return res
        end
    end
    update_global_properties!(Val(:R1), res, glob)
    update_global_properties!(Val(:R2), res, glob)
    update_global_properties!(Val(:off_resonance), res, glob)
    return res
end

function update_global_properties!(::Val{S}, to_update::MRIProperties, glob::GlobalProperties) where {S}
    if isnan(getproperty(to_update, S))
        setproperty!(to_update, S, getproperty(glob, S))
    end
end

function update_stuck_properties!(properties::MRIProperties, group::FixedObstructionGroup, index::Int)
    properties.R1 = get_surface_value(Val(:R1), group, index)
    properties.R2 = get_surface_value(Val(:R2), group, index)
    properties.off_resonance = get_surface_value(Val(:off_resonance), group, index)
end

function get_surface_value(::Val{S}, group::FixedObstructionGroup, index::Int) where {S}
    if ~hasproperty(group.surface, S)
        return NaN
    end
    value = getproperty(group.surface, S)
    if value isa Vector
        return value[index]
    else
        return value
    end
end

function update_inside_properties!(properties::MRIProperties, group::FixedObstructionGroup, position::SVector{3, Float64}, stuck_to::Reflection)
    for index in isinside(group, position, stuck_to)
        update_volume_value!(Val(:R1), properties, group, index)
        update_volume_value!(Val(:R2), properties, group, index)
        update_volume_value!(Val(:off_resonance), properties, group, index)
    end
end

function update_volume_value!(::Val{S}, properties::MRIProperties, group::FixedObstructionGroup, index::Int) where {S}
    if hasproperty(group.volume, S) && isnan(getproperty(properties, S))
        value = getproperty(group.volume, S)
        use_value = value isa Vector ? value[index] : value
        setproperty!(properties, S, use_value)
    end
end

for symbol in (:R1, :R2, :off_resonance)
    @eval begin
        function $symbol(position::AbstractVector, geometry::FixedGeometry, properties::GlobalProperties=GlobalProperties(), stuck_to::Reflection=empty_reflection)
            MRIProperties(geometry, geometry, properties, SVector{3, Float64}(position), stuck_to).$symbol
        end
    end
end


for symbol in (:permeability, :surface_relaxivity, :dwell_time, :surface_density)
    @eval begin
        """
            $($symbol)(properties, )
            $($symbol)(geometry, properties, reflection)
        
        Returns the $($symbol) experienced by the spin hitting the surface represented by a [`Reflection`](@ref).
        The `geometry` has to be a [`FixedGeometry`](@ref).
        """
        function $symbol(group::FixedObstructionGroup, properties::GlobalProperties, has_hit=nothing)
            if hasproperty(group.surface, $(QuoteNode(symbol)))
                local_value = group.surface.$symbol
                if isnothing(has_hit)
                    return local_value
                end
                if local_value isa Vector
                    the_value = local_value[has_hit]
                    if ~isnothing(the_value)
                        return the_value::typeof($symbol(properties))
                    end
                else
                    return local_value
                end
            end
            return $symbol(properties)
        end

        function $symbol(geometry::FixedGeometry, properties::GlobalProperties, has_hit::Tuple{Int, Int})
            $symbol(geometry[has_hit[1]], properties, has_hit[2])
        end

        function $symbol(geometry::FixedGeometry, properties::GlobalProperties, has_hit::Union{Intersection, Reflection})
            $symbol(geometry, properties, (has_hit.geometry_index, has_hit.obstruction_index))
        end
    end
end

"""
    max_timestep_sticking(geometry, properties, diffusivity)

Returns the maximum timestep that can be used while keeping [`stick_probability`](@ref) lower than 1.
"""
function max_timestep_sticking(geometry::FixedGeometry, properties::GlobalProperties, diffusivity::Number)
    minimum([max_timestep_sticking(group, properties, diffusivity) for group in geometry])
end

max_timestep_sticking(geometry::FixedGeometry{0}, properties::GlobalProperties, diffusivity::Number) = Inf

function max_timestep_sticking(group::FixedObstructionGroup, properties::GlobalProperties, diffusivity::Number)
    return maximum(stick_probability.(surface_density(group, properties), dwell_time(group, properties), diffusivity, 1))^-2
end

end