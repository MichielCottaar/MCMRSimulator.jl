module FixProperties

import ...User.Obstructions: Walls, Cylinders, Spheres, Annuli, Mesh, isglobal
import ...Internal.Properties: GeometryLeafProperties, GeometryVectorProperties, GeometryTupleProperties

const _volume_fields = (:R1, :R2, :off_resonance)
const _surface_fields = (:R1, :R2, :off_resonance, :permeability, :dwell_time, :surface_relaxation, :density)

function _property_values(group, field_name::Symbol, category::Symbol, default::Float64)
    user_field_name = field_name === :surface_relaxation ? :relaxation : field_name
    field_value = getproperty(group, Symbol(user_field_name, "_", category))

    if isglobal(field_value)
        value = isnothing(field_value.value) ? default : field_value.value
        return Float64(value)
    end

    values = [
        Float64(isnothing(value) ? default : value)
        for value in field_value.value
    ]
    values
end

function _property_node(values::Float64)
    GeometryLeafProperties(values)
end

function _property_node(values::Vector{Float64})
    GeometryVectorProperties(values)
end

function _property_node(group, field_name::Symbol, category::Symbol, default::Float64)
    _property_node(_property_values(group, field_name, category, default))
end

function _volume_properties(group, category::Symbol)
    values = map(_volume_fields) do field_name
        _property_node(group, field_name, category, 0.)
    end
    NamedTuple{_volume_fields}(Tuple(values))
end

function _surface_properties(
    group,
    category::Symbol;
    permeability=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    density=0.,
)
    defaults = (
        R1 = Float64(0.),
        R2 = Float64(0.),
        off_resonance = Float64(0.),
        permeability = Float64(permeability),
        dwell_time = Float64(dwell_time),
        surface_relaxation = Float64(surface_relaxation),
        density = Float64(density),
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

function _zero_surface(
    ; permeability=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    density=0.,
)
    (
        R1 = GeometryLeafProperties(0.0),
        R2 = GeometryLeafProperties(0.0),
        off_resonance = GeometryLeafProperties(0.0),
        permeability = GeometryLeafProperties(Float64(permeability)),
        dwell_time = GeometryLeafProperties(Float64(dwell_time)),
        surface_relaxation = GeometryLeafProperties(Float64(surface_relaxation)),
        density = GeometryLeafProperties(Float64(density)),
    )
end

function fix_properties(
    group::Union{Cylinders, Spheres, Mesh};
    permeability=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    density=0.,
)
    (
        volume = _volume_properties(group, :inside),
        surface = _surface_properties(
            group,
            :surface;
            permeability,
            dwell_time,
            surface_relaxation,
            density,
        ),
    )
end

function fix_properties(
    group::Walls;
    permeability=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    density=0.,
)
    (
        volume = _zero_volume(),
        surface = _surface_properties(
            group,
            :surface;
            permeability,
            dwell_time,
            surface_relaxation,
            density,
        ),
    )
end

function fix_properties(
    group::Annuli;
    permeability=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    density=0.,
)
    inner_volume_values = NamedTuple{_volume_fields}(
        Tuple(_property_values(group, field, :inner_volume, 0.) for field in _volume_fields),
    )
    outer_volume_values = NamedTuple{_volume_fields}(
        Tuple(_property_values(group, field, :outer_volume, 0.) for field in _volume_fields),
    )
    inner_volume = NamedTuple{_volume_fields}(
        Tuple(
            _property_node(
                getproperty(inner_volume_values, field) .-
                getproperty(outer_volume_values, field),
            )
            for field in _volume_fields
        ),
    )
    outer_volume = NamedTuple{_volume_fields}(
        Tuple(_property_node(getproperty(outer_volume_values, field)) for field in _volume_fields),
    )
    inner_surface = _surface_properties(
        group,
        :inner_surface;
        permeability,
        dwell_time,
        surface_relaxation,
        density,
    )
    outer_surface = _surface_properties(
        group,
        :outer_surface;
        permeability,
        dwell_time,
        surface_relaxation,
        density,
    )

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
