module Fix

include("fix_base_geometry.jl")
include("fix_transformations.jl")
include("fix_properties.jl")
include("fix_susceptibility.jl")

import ..User.Obstructions: ObstructionGroup, Walls, Cylinders, Spheres, Annuli
import ..InternalNew.PhysicalGeometries: PhysicalGeometry
import ..InternalNew.PhysicalGeometries.Groups: GeometryTuple
import ..InternalNew.Properties: GeometryTupleProperties
import ..InternalNew: FixedGeometry

import .FixBaseGeometry: fix_base_geometry
import .FixTransformations: fix_transformations
import .FixProperties: fix_properties
import .FixSusceptibility: fix_susceptibility

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
    susceptibility = fix_susceptibility(group)
    FixedGeometry{typeof(physical_geometry), typeof(properties.volume), typeof(properties.surface), typeof(susceptibility)}(
        physical_geometry,
        properties.volume,
        properties.surface,
        susceptibility,
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
    isempty(geometries) && return FixedGeometry{typeof(GeometryTuple{3}(())), Nothing, Nothing, Tuple{}}(
        GeometryTuple{3}(()),
        nothing,
        nothing,
        (),
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
    susceptibility = fix_susceptibility(geometries)
    FixedGeometry{typeof(physical_geometry), typeof(volume), typeof(surface), typeof(susceptibility)}(
        physical_geometry,
        volume,
        surface,
        susceptibility,
    )
end

end
