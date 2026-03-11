"""
Support to read pulseq files.
"""
module Pulseq

include("types.jl")
include("extensions.jl")
include("basic_parsers.jl")
include("sections_io.jl")
include("components.jl")
include("parsers/parsers.jl")
include("parse_sections.jl")
include("timings.jl")

import .Types: PulseqSequence, PulseqBlock, PulseqSection, PulseqRFPulse, PulseqGradient, AnyPulseqComponent, AnyPulseqGradient, PulseqShape, PulseqExtension, PulseqExtensionDefinition, PulseqADC
import .Extensions: parse_extension, get_extension_name, add_extension_definition!
import .Timings: duration, adc_sample_times, gradient_waveform, rf_pulses


"""
    read_pulseq(IO)

Reads a sequence from a pulseq file (http://pulseq.github.io/).
Pulseq files can be produced using matlab (http://pulseq.github.io/) or python (https://pypulseq.readthedocs.io/en/master/).
"""
function read_pulseq(io::IO)
    sections = SectionsIO.parse_pulseq_sections(io)
    return ParseSections.parse_all_sections(sections)
end

read_pulseq(filename::AbstractString) = open(read_pulseq, filename)

"""
    write_pulseq(IO, sequence)

Writes a sequence to an output IO file.
"""
function write_pulseq(io::IO, sequence::Types.PulseqSequence)
    if sequence.version < v"1.4"
        error("Can only write to pulseq version 1.4 or later.")
    end
    sections = ParseSections.gen_all_sections(sequence)
    for key in [:version, :definitions, :blocks, :rf, :gradients, :trap, :adc, :extensions, :shapes]
        if length(sections[key].content) == 0
            continue
        end
        SectionsIO.write_pulseq_section(io, sections[key])
    end
end

end