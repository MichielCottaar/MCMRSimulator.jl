"""Load geometry from JSON, PLY, or SWC files."""
module LoadGeometry

import ..JSON: read_geometry_json
import ..LoadMesh: load_mesh
import ..LoadSWC: read_swc

function _format(format)
    format = lowercase(String(format))
    startswith(format, ".") && (format = format[2:end])
    Symbol(format)
end

function _read_geometry(io::IO, format; swc_as_spheres=false, kwargs...)
    format = _format(format)
    if format == :json
        isempty(kwargs) || throw(ArgumentError("Keyword arguments are not supported for JSON geometries"))
        return read_geometry_json(io)
    elseif format == :ply
        return load_mesh(io; kwargs...)
    elseif format == :swc
        return read_swc(io; swc_as_spheres=swc_as_spheres, kwargs...)
    end
    throw(ArgumentError("Unsupported geometry format '$format'. Expected :json, :ply, or :swc."))
end

"""
    read_geometry(io::IO; format=:json, kwargs...)

Read geometry from an open stream. Since streams do not have a filename
extension, `format` must be provided for non-JSON input.
"""
function read_geometry(io::IO; format=:json, kwargs...)
    _read_geometry(io, format; kwargs...)
end

"""
    read_geometry(filename::AbstractString; format=nothing, kwargs...)

Read geometry from a JSON, PLY, or SWC file. When `format` is omitted, the
format is inferred from the filename extension.
"""
function read_geometry(filename::AbstractString; format=nothing, kwargs...)
    stripped = strip(filename)
    if isnothing(format) && (startswith(stripped, "{") || startswith(stripped, "["))
        return read_geometry_json(filename)
    end

    selected_format = isnothing(format) ? splitext(filename)[2] : format
    if selected_format isa AbstractString && isempty(selected_format)
        throw(ArgumentError("Could not infer geometry format from '$filename'; provide format explicitly"))
    end
    selected_format = _format(selected_format)
    open(filename, "r") do io
        _read_geometry(io, selected_format; kwargs...)
    end
end

end
