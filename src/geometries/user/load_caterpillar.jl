"""Load overlapping sphere geometries written by CATERPillar."""
module LoadCaterpillar

import StaticArrays: SVector
import ..Obstructions: Spheres

const CATERPILLAR_HEADER = [
    "cell_type", "cell_id", "component", "component_id", "x", "y", "z",
    "inner_radius", "outer_radius",
]

struct CaterpillarRow
    cell_type::String
    cell_id::Int
    component::String
    component_id::Int
    position::SVector{3, Float64}
    inner_radius::Float64
    outer_radius::Float64
end

is_caterpillar_header(line) = lowercase.(split(strip(line))) == CATERPILLAR_HEADER

function _parse_row(line, line_number)
    fields = split(line)
    length(fields) == length(CATERPILLAR_HEADER) || throw(ArgumentError(
        "CATERPillar line $line_number must contain 9 columns, found $(length(fields))",
    ))
    cell_id, component_id = try
        parse(Int, fields[2]), parse(Int, fields[4])
    catch
        throw(ArgumentError("Could not parse CATERPillar IDs on line $line_number"))
    end
    values = try
        parse.(Float64, fields[5:9])
    catch
        throw(ArgumentError("Could not parse CATERPillar coordinates or radii on line $line_number"))
    end
    all(isfinite, values) || throw(ArgumentError(
        "CATERPillar coordinates and radii must be finite on line $line_number",
    ))
    inner_radius, outer_radius = values[4:5]
    inner_radius > 0 && outer_radius > 0 && inner_radius <= outer_radius || throw(ArgumentError(
        "CATERPillar radii must be positive with inner_radius <= outer_radius on line $line_number",
    ))
    CaterpillarRow(
        fields[1], cell_id, fields[3], component_id,
        SVector{3, Float64}(values[1:3]), inner_radius, outer_radius,
    )
end

"""
    read_caterpillar(io_or_filename; kwargs...)

Read a whitespace-delimited CATERPillar output file and return overlapping
[`Spheres`](@ref) groups. `cell_type` values are used only to keep different
cell types in separate groups and are otherwise not restricted.

CATERPillar represents myelinated axon samples with `inner_radius` and
`outer_radius`. If any row belonging to a cell has different radii, all rows
for that cell are included in both an inner-sphere group and an outer-sphere
group. For rows with equal radii in such a cell, `sqrt(eps(Float64))` is added
to the outer radius so the two spheres do not coincide exactly. Magnetic
susceptibility of the imported myelin is not currently supported; the
generated sphere groups use the supplied sphere properties only.

Additional keyword arguments are passed to every generated [`Spheres`](@ref)
group. CATERPillar spheres are always loaded with `overlapping=true`.
"""
function read_caterpillar(io::IO; kwargs...)
    rows = [
        (line_number, strip(line))
        for (line_number, line) in enumerate(eachline(io))
        if !isempty(strip(line)) && !startswith(strip(line), "#")
    ]
    isempty(rows) && throw(ArgumentError("CATERPillar file is empty"))
    is_caterpillar_header(rows[1][2]) || throw(ArgumentError(
        "CATERPillar header must be 'cell_type cell_id component component_id X Y Z inner_radius outer_radius'",
    ))
    cells = Dict{Tuple{String, Int}, Vector{CaterpillarRow}}()
    cell_order = Tuple{String, Int}[]
    for (line_number, line) in rows[2:end]
        row = _parse_row(line, line_number)
        key = (row.cell_type, row.cell_id)
        if !haskey(cells, key)
            cells[key] = CaterpillarRow[]
            push!(cell_order, key)
        end
        push!(cells[key], row)
    end
    isempty(cells) && throw(ArgumentError("CATERPillar file contains no data rows"))

    positions = Dict{Tuple{String, Symbol}, Vector{SVector{3, Float64}}}()
    radii = Dict{Tuple{String, Symbol}, Vector{Float64}}()
    layer_order = Tuple{String, Symbol}[]
    function add_layer!(key, row, radius)
        if !haskey(positions, key)
            positions[key] = SVector{3, Float64}[]
            radii[key] = Float64[]
            push!(layer_order, key)
        end
        push!(positions[key], row.position)
        push!(radii[key], radius)
    end

    for cell_key in cell_order
        cell_type = first(cell_key)
        cell_rows = cells[cell_key]
        has_shell = any(row.inner_radius != row.outer_radius for row in cell_rows)
        if has_shell
            for row in cell_rows
                add_layer!((cell_type, :inner), row, row.inner_radius)
                outer_radius = row.outer_radius
                row.inner_radius == row.outer_radius && (outer_radius += sqrt(eps(Float64)))
                add_layer!((cell_type, :outer), row, outer_radius)
            end
        else
            for row in cell_rows
                add_layer!((cell_type, :ordinary), row, row.outer_radius)
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
