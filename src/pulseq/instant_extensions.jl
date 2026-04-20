"""
This module contains the code for handling the InstantPulse and InstantGradient extensions in Pulseq files. 
It defines how to parse these extensions from a Pulseq file.
"""
module InstantExtensions
import ..Extensions: parse_extension, get_extension_name, add_extension_definition!, PulseqExtensionDefinition
import StaticArrays: SVector


struct InstantPulse
    flip_angle::Float64
    phase::Float64
end

# I/O of InstantPulse
function parse_extension(ext::PulseqExtensionDefinition{:InstantPulse})
    mapping = Dict{Int, Tuple{Float64, InstantPulse}}()
    for line in ext.content
        (id, delay, flip_angle, phase) = parse.((Int, Float64, Float64, Float64), split(line))
        mapping[id] = (delay, InstantPulse(flip_angle, phase))
    end
    return mapping
end

get_extension_name(::Tuple{<:Number, InstantPulse}) = :InstantPulse

function add_extension_definition!(content::Vector{String}, obj::Tuple{Number, InstantPulse})
    to_store = (obj[1], obj[2].flip_angle, obj[2].phase)
    for line in content
        (id, this_line...) = parse.((Int, Float64, Float64, Float64), split(line))
        if all(to_store .≈ this_line)
            return id
        end
    end
    push!(content, "$(length(content) + 1) " * join(string.(to_store), " "))
    return length(content)
end

struct InstantGradient
    qvec::SVector{3, Float64}
end

# I/O of InstantGradient
function parse_extension(ext::PulseqExtensionDefinition{:InstantGradient})
    mapping = Dict{Int, Tuple{Float64, InstantGradient}}()
    for line in ext.content
        (id, delay, qvec...) = parse.((Int, Float64, Float64, Float64, Float64), split(line))
        mapping[id] = (delay, InstantGradient(SVector{3, Float64}(qvec...)))
    end
end

get_extension_name(::Tuple{<:Number, InstantGradient}) = :InstantGradient

function add_extension_definition!(content::Vector{String}, obj::Tuple{<:Number, <:InstantGradient})
    to_store = (obj[1], variables.qvec(obj[2])...)
    for line in content
        (id, this_line...) = parse.((Int, Float64, Float64, Float64, Float64), split(line))
        if all(to_store .≈ this_line)
            return id
        end
    end
    push!(content, "$(length(content) + 1) " * join(string.(to_store), " "))
    return length(content)
end

end