module Fix

include("fix_base_geometry.jl")
include("fix_transformations.jl")
include("fix_properties.jl")
include("fix_susceptibility.jl")

import ..User.Obstructions: ObstructionGroup, IndexedObstruction, Walls, Cylinders, Spheres, Annuli, BendyCylinder, Mesh
import ..Internal.PhysicalGeometries: PhysicalGeometry
import ..Internal: SizeScaleOverride
import ..Internal.PhysicalGeometries.Groups: GeometryTuple
import ..Internal.Properties: GeometryTupleProperties
import ..Internal: FixedGeometry

import .FixBaseGeometry: fix_base_geometry, annuli_size_scale
import .FixTransformations: fix_transformations
import .FixProperties: fix_properties, _zero_volume, _zero_surface
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

function fix(group::BendyCylinder; kwargs...)
    fix(Mesh(group); kwargs...)
end

function fix(obstruction::IndexedObstruction; kwargs...)
    group = obstruction.group
    values = (; (key => getproperty(obstruction, key) for key in propertynames(group))...)
    fix(ObstructionGroup(group.type; number=1, values...); kwargs...)
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
    size_scale = group.size_scale.value
    if isnothing(size_scale) && group isa Annuli
        size_scale = annuli_size_scale(group)
    end
    if !isnothing(size_scale)
        physical_geometry = SizeScaleOverride(physical_geometry, size_scale)
    end
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
    if isempty(geometries)
        volume = _zero_volume()
        surface = _zero_surface(; permeability, dwell_time, surface_relaxation, density)
        return FixedGeometry{typeof(GeometryTuple{3}(())), typeof(volume), typeof(surface), Tuple{}}(
            GeometryTuple{3}(()),
            volume,
            surface,
            (),
        )
    end

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

fix(geometry::Tuple; kwargs...) = fix(collect(geometry); kwargs...)

end
