"""
Support to read pulseq files.
"""
module Pulseq

include("types.jl")
include("extensions.jl")
include("instant_extensions.jl")
include("basic_parsers.jl")
include("sections_io.jl")
include("components.jl")
include("parsers/parsers.jl")
include("parse_sections.jl")
include("timings.jl")

import .Types: PulseqSequence, PulseqBlock, PulseqSection, PulseqRFPulse, PulseqGradient, AnyPulseqComponent, AnyPulseqGradient, PulseqShape, PulseqExtension, PulseqExtensionDefinition, PulseqADC
import .Extensions: parse_extension, get_extension_name, add_extension_definition!
import .InstantExtensions: InstantPulse, InstantGradient
import .Timings: duration, adc_sample_times, gradient_waveform, rf_pulses


"""
    read_pulseq(IO)

Reads a sequence from a pulseq file (http://pulseq.github.io/).
Pulseq files can be produced using matlab (http://pulseq.github.io/) or python (https://pypulseq.readthedocs.io/en/master/).
"""
function read_pulseq(io::IO; TR=nothing)
    sections = SectionsIO.parse_pulseq_sections(io)
    sequence = ParseSections.parse_all_sections(sections)
    return _set_repetition_time(sequence, TR)
end

read_pulseq(filename::AbstractString; TR=nothing) = open(filename) do io
    read_pulseq(io; TR=TR)
end

function _set_repetition_time(sequence::Types.PulseqSequence, TR)
    if isnothing(TR)
        return sequence
    end
    TR_seconds = Float64(TR) / 1000
    sequence_duration = duration(sequence, :second)
    if sequence_duration > TR_seconds
        error("The sequence duration ($sequence_duration s) is longer than the specified TR ($TR_seconds s).")
    end
    return Types.PulseqSequence(sequence.version, sequence.definitions, sequence.blocks, TR_seconds)
end

"""
    write_pulseq(IO, sequence)

Writes a sequence to an output IO file.
"""
function write_pulseq(io::IO, sequence::Types.PulseqSequence)
    if sequence.version < v"1.5"
        error("Can only write to pulseq version 1.5 or later.")
    end
    sections = ParseSections.gen_all_sections(sequence)
    for key in [:version, :definitions, :blocks, :rf, :gradients, :trap, :adc, :extensions, :shapes]
        if length(sections[key].content) == 0
            continue
        end
        SectionsIO.write_pulseq_section(io, sections[key])
    end
end

write_pulseq(filename::AbstractString, sequence::Types.PulseqSequence) = open(filename, "w") do io
    write_pulseq(io, sequence)
end

end
