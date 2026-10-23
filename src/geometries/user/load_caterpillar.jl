"""Load overlapping sphere geometries written by CATERPillar."""
module LoadCaterpillar

import StaticArrays: SVector
import CSV
import ..Obstructions: Spheres

const CATERPILLAR_HEADER = [
    "cell_type", "cell_id", "component", "component_id", "x", "y", "z",
    "inner_radius", "outer_radius",
]
const CATERPILLAR_HEADER_WITH_PARENT = [
    "cell_type", "cell_id", "component", "component_id", "parent_component_id",
    "x", "y", "z", "inner_radius", "outer_radius",
]

function is_caterpillar_header(line)
    header = lowercase.(split(strip(line)))
    header == CATERPILLAR_HEADER || header == CATERPILLAR_HEADER_WITH_PARENT
end

function _is_valid_header(names)
    names = String.(names)
    normalized = lowercase.(names)
    if normalized == CATERPILLAR_HEADER
        return names[5:7] == ["X", "Y", "Z"]
    elseif normalized == CATERPILLAR_HEADER_WITH_PARENT
        return names[6:8] == ["X", "Y", "Z"]
    end
    false
end

"""
    read_caterpillar(io_or_filename; kwargs...)

Read a whitespace-delimited CATERPillar output file and return overlapping
[`Spheres`](@ref) groups. Each `(cell_type, cell_id)` produces its own sphere
group, or separate inner and outer groups when the cell has a shell. `cell_type`
values are otherwise not restricted.

CATERPillar represents myelinated axon samples with `inner_radius` and
`outer_radius`. If any row belonging to a cell has different radii, all rows
for that cell are included in both an inner-sphere group and an outer-sphere
group. For rows with equal radii in such a cell, `sqrt(eps(Float64))` is added
to the outer radius so the two spheres do not coincide exactly. Magnetic
susceptibility of the imported myelin is not currently supported; the
generated sphere groups use the supplied sphere properties only.

Both the nine-column format and the ten-column format with an integer
`parent_component_id` between `component_id` and `X` are supported.
The parent ID is validated but does not alter the generated sphere geometry.

Additional keyword arguments are passed to every generated [`Spheres`](@ref)
group. CATERPillar spheres are always loaded with `overlapping=true`.
"""
function read_caterpillar(io::IO; kwargs...)
    table = try
        CSV.File(io; comment="#", on_error=:collect)
    catch err
        throw(ArgumentError("Could not parse CATERPillar table: $(sprint(showerror, err))"))
    end
    read_caterpillar(table; kwargs...)
end

function read_caterpillar(table::CSV.File; kwargs...)
    _is_valid_header(propertynames(table)) || throw(ArgumentError(
        "CATERPillar header must contain the standard columns, optionally including parent_component_id",
    ))
    has_parent_component = :parent_component_id in propertynames(table)
    integer_columns = has_parent_component ?
        (:cell_id, :component_id, :parent_component_id) : (:cell_id, :component_id)
    all(name -> Base.nonmissingtype(eltype(getproperty(table, name))) <: Integer, integer_columns) ||
        throw(ArgumentError("CATERPillar IDs must be integers"))
    all(name -> Base.nonmissingtype(eltype(getproperty(table, name))) <: Real,
        (:X, :Y, :Z, :inner_radius, :outer_radius)) ||
        throw(ArgumentError("CATERPillar coordinates and radii must be numeric"))
    read_caterpillar(
        table.cell_type, table.cell_id, table.X, table.Y, table.Z,
        table.inner_radius, table.outer_radius; kwargs...,
    )
end

function read_caterpillar(cell_type, cell_id, X, Y, Z, inner_radius, outer_radius; kwargs...)
    cells = Dict{Tuple{String, Int}, Vector{Tuple{SVector{3, Float64}, Float64, Float64}}}()
    cell_order = Tuple{String, Int}[]
    for i in eachindex(X)
        line_number = i + 1
        fields = (cell_type[i], cell_id[i], X[i], Y[i], Z[i], inner_radius[i], outer_radius[i])
        any(ismissing, fields) && throw(ArgumentError(
            "Missing or invalid CATERPillar field on line $line_number",
        ))
        all(isfinite, (X[i], Y[i], Z[i], inner_radius[i], outer_radius[i])) || throw(ArgumentError(
            "CATERPillar coordinates and radii must be finite on line $line_number",
        ))
        inner_radius[i] > 0 && outer_radius[i] > 0 && inner_radius[i] <= outer_radius[i] || throw(ArgumentError(
            "CATERPillar radii must be positive with inner_radius <= outer_radius on line $line_number",
        ))
        position = SVector{3, Float64}(X[i], Y[i], Z[i])
        key = (String(cell_type[i]), Int(cell_id[i]))
        if !haskey(cells, key)
            cells[key] = Tuple{SVector{3, Float64}, Float64, Float64}[]
            push!(cell_order, key)
        end
        push!(cells[key], (position, Float64(inner_radius[i]), Float64(outer_radius[i])))
    end
    isempty(cells) && throw(ArgumentError("CATERPillar file contains no data rows"))

    positions = Dict{Tuple{String, Int, Symbol}, Vector{SVector{3, Float64}}}()
    radii = Dict{Tuple{String, Int, Symbol}, Vector{Float64}}()
    layer_order = Tuple{String, Int, Symbol}[]
    function add_layer!(key, position, radius)
        if !haskey(positions, key)
            positions[key] = SVector{3, Float64}[]
            radii[key] = Float64[]
            push!(layer_order, key)
        end
        push!(positions[key], position)
        push!(radii[key], radius)
    end

    for cell_key in cell_order
        cell_type, cell_id = cell_key
        cell_rows = cells[cell_key]
        has_shell = any(inner_radius != outer_radius for (_, inner_radius, outer_radius) in cell_rows)
        if has_shell
            for (position, inner_radius, outer_radius) in cell_rows
                add_layer!((cell_type, cell_id, :inner), position, inner_radius)
                outer_radius == inner_radius && (outer_radius += sqrt(eps(Float64)))
                add_layer!((cell_type, cell_id, :outer), position, outer_radius)
            end
        else
            for (position, _, outer_radius) in cell_rows
                add_layer!((cell_type, cell_id, :ordinary), position, outer_radius)
            end
        end
    end

    haskey(kwargs, :overlapping) && throw(ArgumentError(
        "overlapping is fixed to true for CATERPillar geometries",
    ))
    [
        Spheres(; position=positions[key], radius=radii[key], overlapping=true, kwargs...)
        for key in layer_order
    ]
end

read_caterpillar(filename::AbstractString; kwargs...) =
    open(io -> read_caterpillar(io; kwargs...), filename; read=true)

end
