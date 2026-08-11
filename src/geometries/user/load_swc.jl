"""
    SWCNode
    SWCFile
    read_swc

Representation and loading of files in the standard SWC neuron morphology
format. Coordinates and radii retain the micrometre units used by SWC files.
"""
module LoadSWC

import StaticArrays: SVector
import ..Obstructions: Spheres

"""One row of an SWC file."""
struct SWCNode
    id::Int
    type_id::Int
    position::SVector{3, Float64}
    radius::Float64
    parent_id::Int
end

"""
An SWC file represented by its comments and ordered node table.

Comments before the first node are stored in `header`; comments after the
first node are stored in `footer`.
"""
struct SWCFile
    header::Vector{String}
    nodes::Vector{SWCNode}
    footer::Vector{String}
end

function _parse_node(line, line_number)
    fields = split(strip(line))
    if length(fields) != 7
        throw(ArgumentError("SWC line $line_number must contain 7 columns, found $(length(fields))"))
    end

    try
        return SWCNode(
            parse(Int, fields[1]),
            parse(Int, fields[2]),
            SVector{3, Float64}(parse.(Float64, fields[3:5])),
            parse(Float64, fields[6]),
            parse(Int, fields[7]),
        )
    catch error
        error isa ArgumentError || rethrow()
        throw(ArgumentError("Could not parse SWC line $line_number: $line"))
    end
end

function _validate_node(node, line_number, ids)
    node.id > 0 || throw(ArgumentError("SWC node ID on line $line_number must be positive"))
    node.type_id >= 0 || throw(ArgumentError("SWC type ID on line $line_number must be non-negative"))
    node.radius >= 0 || throw(ArgumentError("SWC radius on line $line_number must be non-negative"))
    node.id ∉ ids || throw(ArgumentError("Duplicate SWC node ID $(node.id) on line $line_number"))
    if node.parent_id != -1 && node.parent_id ∉ ids
        throw(ArgumentError("SWC parent ID $(node.parent_id) on line $line_number is not an earlier node"))
    end
    push!(ids, node.id)
end

"""
    read_swc_raw(filename::AbstractString)
    read_swc_raw(io::IO)

Read an SWC file from `io` or `filename`.

Blank lines are ignored and comment lines are retained as header or footer lines.
"""
function read_swc_raw(io::IO)
    header = String[]
    footer = String[]
    nodes = SWCNode[]
    ids = Set{Int}()

    for (line_number, line) in enumerate(eachline(io))
        stripped = strip(line)
        isempty(stripped) && continue
        if startswith(stripped, '#')
            push!(isempty(nodes) ? header : footer, line)
            continue
        end

        node = _parse_node(line, line_number)
        if isempty(nodes) && node.parent_id != -1
            throw(ArgumentError("The first SWC node must have parent ID -1"))
        end
        _validate_node(node, line_number, ids)
        push!(nodes, node)
    end

    isempty(nodes) && throw(ArgumentError("SWC file contains no nodes"))
    count(node -> node.parent_id == -1, nodes) == 1 ||
        throw(ArgumentError("SWC file must contain exactly one root node"))
    SWCFile(header, nodes, footer)
end

function read_swc_raw(filename::AbstractString)
    open(read_swc_raw, filename; read=true)
end


"""
    read_swc(swc_file; swc_as_spheres=false, kwargs...)

Read an SWC file and return a `Spheres` object with the node positions and radii.

See [`Spheres`](@ref) for the available keyword arguments. 

In the future, this function may be extended to return a more complete representation of the SWC file, including the links between spheres. 
For now, it only returns the spheres themselves with the `overlapping` keyword set to `true` by default.

For forwards compatibility, the `swc_as_spheres` keyword argument is provided, but it must be set to `true` to avoid an error.
In the future, this argument may be removed when the function is extended to return a more complete representation of the SWC file.
"""
function read_swc(swc_file::SWCFile; swc_as_spheres=false, kwargs...)
    if !swc_as_spheres
        throw(ArgumentError("Reading SWC files as geometries including the links between spheres is not supported yet. Set `swc_as_spheres` to true if you want to load them as individual overlapping spheres without linking cylinders."))
    end
    radii = [node.radius for node in swc_file.nodes]
    positions = [node.position for node in swc_file.nodes]
    return Spheres(; position=positions, radius=radii, overlapping=swc_as_spheres, kwargs...)
end

read_swc(in_file::Union{IO, AbstractString}; kwargs...) = read_swc(read_swc_raw(in_file); kwargs...)

end
