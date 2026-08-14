module FixProperties

import ...User.Obstructions: Walls, Cylinders, Spheres, Annuli, isglobal
import ...InternalNew.Properties: GeometryLeafProperties, GeometryVectorProperties, GeometryTupleProperties

const _volume_fields = (:R1, :R2, :off_resonance)
const _surface_fields = (:R1, :R2, :off_resonance, :permeability, :dwell_time, :surface_relaxation, :density)

function _property_node(group, field_name::Symbol, category::Symbol, default::Float64)
    user_field_name = field_name === :surface_relaxation ? :relaxation : field_name
    field_value = getproperty(group, Symbol(user_field_name, "_", category))

    if isglobal(field_value)
        value = isnothing(field_value.value) ? default : field_value.value
        return GeometryLeafProperties(Float64(value))
    end

    values = [
        Float64(isnothing(value) ? default : value)
        for value in field_value.value
    ]
    GeometryVectorProperties(values)
end

function _volume_properties(group, category::Symbol)
    values = map(_volume_fields) do field_name
        _property_node(group, field_name, category, 0.)
    end
    NamedTuple{_volume_fields}(Tuple(values))
end

function _surface_properties(group, category::Symbol; kwargs...)
    defaults = (
        R1 = Float64(0.),
        R2 = Float64(0.),
        off_resonance = Float64(0.),
        permeability = Float64(get(kwargs, :permeability, 0.0)),
        dwell_time = Float64(get(kwargs, :dwell_time, 0.0)),
        surface_relaxation = Float64(get(kwargs, :relaxation, 0.0)),
        density = Float64(get(kwargs, :density, 0.0)),
    )
    values = map(_surface_fields) do field_name
        _property_node(group, field_name, category, getproperty(defaults, field_name))
    end
    NamedTuple{_surface_fields}(Tuple(values))
end

function _zero_volume()
    (
        R1 = GeometryLeafProperties(0.0),
        R2 = GeometryLeafProperties(0.0),
        off_resonance = GeometryLeafProperties(0.0),
    )
end

function fix_properties(group::Union{Cylinders, Spheres}; kwargs...)
    (
        volume = _volume_properties(group, :inside),
        surface = _surface_properties(group, :surface; kwargs...),
    )
end

function fix_properties(group::Walls; kwargs...)
    (
        volume = _zero_volume(),
        surface = _surface_properties(group, :surface; kwargs...),
    )
end

function fix_properties(group::Annuli; kwargs...)
    inner_volume = _volume_properties(group, :inner_volume)
    outer_volume = _volume_properties(group, :outer_volume)
    inner_surface = _surface_properties(group, :inner_surface; kwargs...)
    outer_surface = _surface_properties(group, :outer_surface; kwargs...)

    volume = NamedTuple{_volume_fields}(
        Tuple(GeometryTupleProperties((getproperty(inner_volume, field), getproperty(outer_volume, field)))
              for field in _volume_fields),
    )
    surface = NamedTuple{_surface_fields}(
        Tuple(GeometryTupleProperties((getproperty(inner_surface, field), getproperty(outer_surface, field)))
              for field in _surface_fields),
    )
    return (
        volume = volume,
        surface = surface,
    )
end

end
