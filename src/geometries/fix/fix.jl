module Fix

include("fix_base_geometry.jl")
include("fix_transformations.jl")
include("fix_properties.jl")

import ..User.Obstructions: ObstructionGroup, Walls, Cylinders, Spheres, Annuli
import ..InternalNew.PhysicalGeometries: PhysicalGeometry
import ..InternalNew.PhysicalGeometries.Groups: GeometryTuple
import ..InternalNew.Properties: GeometryTupleProperties
import ..InternalNew: FixedGeometry

import .FixBaseGeometry: fix_base_geometry
import .FixTransformations: fix_transformations
import .FixProperties: fix_properties

function fix(
    geometry::FixedGeometry;
    permeability=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    density=0.,
)
    (permeability == 0. && dwell_time == 0. && surface_relaxation == 0. && density == 0.) ||
        throw(ArgumentError("cannot apply fix keywords to an already fixed geometry"))
    geometry
end

function fix(
    group::ObstructionGroup;
    permeability=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    density=0.,
)
    base_geometry = fix_base_geometry(group)
    physical_geometry = fix_transformations(group, base_geometry)
    properties = fix_properties(
        group;
        permeability,
        dwell_time,
        surface_relaxation,
        density,
    )
    FixedGeometry{typeof(physical_geometry), typeof(properties.volume), typeof(properties.surface)}(
        physical_geometry,
        properties.volume,
        properties.surface,
    )
end

function _merge_properties(fixed, field_name::Symbol)
    first_properties = getproperty(first(fixed), field_name)
    fields = propertynames(first_properties)
    NamedTuple{fields}(
        Tuple(
            GeometryTupleProperties(Tuple(getproperty(getproperty(geometry, field_name), field) for geometry in fixed))
            for field in fields
        ),
    )
end

function fix(
    geometries::AbstractVector;
    permeability=0.,
    dwell_time=0.,
    surface_relaxation=0.,
    density=0.,
)
    isempty(geometries) && return FixedGeometry{typeof(GeometryTuple{3}(())), Nothing, Nothing}(
        GeometryTuple{3}(()),
        nothing,
        nothing,
    )

    fixed = [
        fix(
            geometry;
            permeability,
            dwell_time,
            surface_relaxation,
            density,
        )
        for geometry in geometries
    ]
    physical_geometry = GeometryTuple(Tuple(geometry.geometry for geometry in fixed))
    volume = _merge_properties(fixed, :volume)
    surface = _merge_properties(fixed, :surface)
    FixedGeometry{typeof(physical_geometry), typeof(volume), typeof(surface)}(
        physical_geometry,
        volume,
        surface,
    )
end

end
