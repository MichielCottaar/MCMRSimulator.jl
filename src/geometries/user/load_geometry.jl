"""Load geometry from JSON, PLY, SWC, or liminal files."""
module LoadGeometry

import ..JSON: read_geometry_json
import ..LoadMesh: load_mesh
import ..LoadSWC: read_swc
import ..LiminalGeometries: LiminalGeometry

function _format(format)
    format = lowercase(String(format))
    startswith(format, ".") && (format = format[2:end])
    Symbol(format)
end

function _first_content_line(io::IO)
    for line in eachline(io)
        stripped = strip(line)
        isempty(stripped) && continue
        startswith(stripped, "#") && continue
        return stripped
    end
    nothing
end

function _detect_format(io::IO)
    mark(io)
    try
        firstline = _first_content_line(io)
        isnothing(firstline) && throw(ArgumentError("Could not detect geometry format from an empty file"))
        first(firstline) in ('{', '[') && return :json
        firstline == "ply" && return :ply
        lowercase(first(split(firstline))) == "liminal" && return :liminal
        :swc
    finally
        reset(io)
    end
end

function _read_liminal(io::IO; base_dir=pwd(), swc_as_spheres=false)
    rows = [(line_number, strip(line)) for (line_number, line) in enumerate(eachline(io)) if
        !isempty(strip(line)) && !startswith(strip(line), "#")]
    isempty(rows) && throw(ArgumentError("Liminal geometry file is empty"))

    header_line, header = first(rows)
    header_fields = split(header)
    length(header_fields) == 2 && lowercase(header_fields[1]) == "liminal" ||
        throw(ArgumentError("Liminal header on line $header_line must be 'liminal <extracellular fraction>'"))
    extracellular_fraction = try
        parse(Float64, header_fields[2])
    catch
        throw(ArgumentError("Could not parse extracellular fraction on line $header_line"))
    end

    geometries = Tuple{Float64, Any}[]
    for (line_number, line) in rows[2:end]
        fields = split(line)
        length(fields) == 2 || throw(ArgumentError(
            "Liminal cell definition on line $line_number must contain '<fraction> <filename>'",
        ))
        fraction = try
            parse(Float64, fields[1])
        catch
            throw(ArgumentError("Could not parse cell fraction on line $line_number"))
        end
        filename = fields[2]
        child_filename = isabspath(filename) ? filename : joinpath(base_dir, filename)
        isfile(child_filename) || throw(ArgumentError(
            "Liminal child geometry file does not exist on line $line_number: $filename",
        ))
        push!(geometries, (fraction, read_geometry(child_filename; swc_as_spheres)))
    end

    LiminalGeometry(; geometries, extracellular_fraction)
end

function _read_geometry(io::IO, format; base_dir=pwd(), swc_as_spheres=false, kwargs...)
    isnothing(format) && (format = _detect_format(io))
    format = _format(format)
    if format == :json
        isempty(kwargs) || throw(ArgumentError("Keyword arguments are not supported for JSON geometries"))
        return read_geometry_json(io)
    elseif format == :ply
        return load_mesh(io; kwargs...)
    elseif format == :swc
        return read_swc(io; swc_as_spheres=swc_as_spheres, kwargs...)
    elseif format == :liminal
        isempty(kwargs) || throw(ArgumentError("Keyword arguments are not supported for liminal geometries"))
        return _read_liminal(io; base_dir, swc_as_spheres)
    end
    throw(ArgumentError("Unsupported geometry format '$format'. Expected :json, :ply, :swc, or :liminal."))
end

"""
    read_geometry(io::IO; format=nothing, kwargs...)

Read geometry from an open stream. If `format` is omitted, the format is
detected from the first non-empty content line. JSON, PLY, SWC, and liminal
geometry files are supported. Use `format` to override content detection.
"""
function read_geometry(io::IO; format=nothing, kwargs...)
    _read_geometry(io, format; kwargs...)
end

"""
    read_geometry(filename::AbstractString; format=nothing, kwargs...)

Read geometry from a JSON, PLY, SWC, or liminal geometry file. When `format`
is omitted, the format is detected from the file contents rather than its
filename extension. Liminal files contain an extracellular volume fraction
and references to child geometry files; relative child paths are resolved
relative to the liminal file.
"""
function read_geometry(filename::AbstractString; format=nothing, kwargs...)
    stripped = strip(filename)
    if isnothing(format) && (startswith(stripped, "{") || startswith(stripped, "["))
        return read_geometry_json(filename)
    end

    open(filename, "r") do io
        _read_geometry(io, format; base_dir=dirname(abspath(filename)), kwargs...)
    end
end

end
