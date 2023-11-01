module Obstructions
include("fields.jl")
include("obstruction_types.jl")
include("obstruction_groups.jl")
import StaticArrays: MVector
import .ObstructionTypes: ObstructionType, fields
import .Fields: Field, FieldValue, isglobal, description
import .ObstructionGroups: ObstructionGroup, IndexedObstruction, key_value_pairs
import ....Methods: get_rotation

function field_to_docs(key::Symbol, field_value::FieldValue{T}) where {T}
    key_txt = replace(String(key), "_"=>"\\_")
    return "- $(key_txt): $(field_to_docs(field_value))\n"
end

function field_to_docs(field_value::FieldValue{T}) where {T}
    field = field_value.field
    as_string = description(field_value)
    if ~isnothing(field.default_value)
        as_string = as_string * "  default value: $(field.default_value)"
    end
    return as_string
end

for obstruction_type in (
    ObstructionType(:Wall; ndim=1, volumes=[]),
    ObstructionType(
        :Cylinder; ndim=2, fields=[
            Field{Float64}(:radius, "Radius of the cylinder.", required=true), 
            Field{Float64}(:g_ratio, "Inner/outer radius used for susceptibility calculation", 1.),
            Field{Float64}(:susceptibility_iso, "Isotropic component of the susceptibility (in ppm).", -0.1),
            Field{Float64}(:susceptibility_aniso, "Ansotropic component of the susceptibility (in ppm).", -0.1),
            Field{Float64}(:lorentz_radius, "Only compute field explicitly for cylinders with this Lorentz radius.", 5.),
        ]),
    ObstructionType(
        :Annulus; plural=:Annuli, ndim=2, surfaces=[:inner_surface, :outer_surface], volumes=[:inner_volume, :outer_volume], fields=[
            Field{Float64}(:inner, "Radius of the inner cylinder.", required=true), 
            Field{Float64}(:outer, "Radius of the outer cylinder.", required=true), 
            Field{Bool}(:myelin, "Whether the annulus is myelinated.", false, required=true), 
            Field{Float64}(:susceptibility_iso, "Isotropic component of the myelin susceptibility (in ppm).", -0.1),
            Field{Float64}(:susceptibility_aniso, "Ansotropic component of the myelin susceptibility (in ppm).", -0.1),
            Field{Float64}(:lorentz_radius, "Only compute field explicitly for a annuli with this Lorentz radius.", 5.),
        ]),
    ObstructionType(
        :Sphere; ndim=3, fields=[
            Field{Float64}(:radius, "Radius of the cylinder.", required=true), 
        ]),
    ObstructionType(
        :Triangle; plural=:Mesh, ndim=3, include_shift=false, fields=[
            Field{MVector{3, Int}}(:triangles, "Each triangle is defined by 3 vertices into the mesh.", required=true),
            Field{Vector{MVector{3, Float64}}}(:vertices, "Positions of the corners of the triangular mesh.", required=true, only_group=true),
            Field{Bool}(:myelin, "Whether the mesh is myelinated.", false, required=true), 
            Field{Float64}(:susceptibility_iso, "Isotropic component of the myelin susceptibility (in ppm).", -0.1),
            Field{Float64}(:susceptibility_aniso, "Ansotropic component of the myelin susceptibility (in ppm).", -0.1),
            Field{Bool}(:save_memory, "If true the internal triangle representation will contain the indices of the vertices rather than the actual corner coordinates. This will save memory, but slow down the algorithm as more memory lookups will have to take place.", true, required=true, only_group=true),
            Field{Float64}(:lorentz_radius, "Only compute field explicitly for triangles with this Lorentz radius.", 5.),
        ]),
    ObstructionType(
        :Ring; plural=:BendyCylinder, ndim=3, include_shift=false, fields=[
            Field{MVector{3, Float64}}(:control_point, "Control points defining the path of the cylinder.", required=true), 
            Field{Float64}(:radius, "Radius at each control point.", required=true), 
            Field{Int}(:nsamples, "Number of mesh vertices along each ring.", 100, only_group=true),
            Field{Bool}(:myelin, "Whether the rings are mylinated.", false, required=true, only_group=true), 
            Field{Float64}(:susceptibility_iso, "Isotropic component of the myelin susceptibility (in ppm).", -0.1),
            Field{Float64}(:susceptibility_aniso, "Ansotropic component of the myelin susceptibility (in ppm).", -0.1),
        ]),
)
    name_string = String(obstruction_type.singular)
    (unique_keys, field_values) = key_value_pairs(obstruction_type, 0)
    unique_field_values = [field_values[key] for key in unique_keys]
    sort!(unique_field_values, by=fv->~fv.field.required)
    fields_string = join(field_to_docs.(unique_keys, unique_field_values))
    @eval const $(obstruction_type.singular) = IndexedObstruction{$(QuoteNode(obstruction_type.plural))}
    @eval const $(obstruction_type.plural) = ObstructionGroup{$(QuoteNode(obstruction_type.plural))}
    constructor = obstruction_type.plural
    @eval begin
        $(constructor)(; kwargs...) = ObstructionGroup($(obstruction_type); kwargs...)
        @doc """
            $($constructor))(; fields...)

        Creates a set of """ * $(name_string) * """ objects.
        Fields can be set using keyword arguments.
        The following fields are available:
        """ * $(fields_string) $constructor
    end
end

end