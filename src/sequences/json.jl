"""
I/O for sequences using JSON format

Functions:
- `write_sequence`
- `read_sequence`
"""
module JSON

import JSON as JSON_pkg
import StaticArrays: SVector
import ..Main: Sequence
import ...Scanners: Scanner
import ..Instants: InstantRFPulse, InstantGradient, Readout
import ..RadioFrequency: RFPulse
import ..Gradients: MRGradients
import ..Shapes: Shape
import ..PulseQ: read_pulseq

"""
    write_sequence([io, ]sequence)

Writes the sequence as a JSON file.
If no IO is provided, the sequence is written to stdout.
"""
write_sequence(object) = write_sequence(stdout, object)

function write_sequence(io::IO, sequence::Sequence) 
    as_string = JSON_pkg.json(sequence)
    print(io, as_string)
end

function write_sequence(filename::AbstractString, sequence)
    open(io -> write_sequence(io, sequence), filename; write=true)
end

"""
    read_sequence(filename/io)

Read sequence from a JSON file.
"""
function read_sequence(io::IO)
    obj = JSON_pkg.parse(io)
    parse_sequence(obj)
end

function read_sequence(filename::AbstractString)
    if startswith(filename, "{")
        obj = JSON_pkg.parse(filename)
        return parse_sequence(obj)
    end
    if endswith(filename, ".seq")
        return read_pulseq(filename)
    end
    open(read_sequence, filename; read=true)
end


function parse_sequence(d::Dict)
    Sequence(;
        TR=d["TR"],
        scanner=parse_scanner(d["scanner"]),
        components=[
            parse_instant.(d["instants"])...,
            parse_pulse.(d["pulses"])...,
            Readout.(d["readout_times"])...,
            parse_gradients.(d["gradients"])...,
        ],
    )
end

function parse_instant(d::Dict)
    if haskey(d, "flip_angle")
        InstantRFPulse(; time=d["time"], flip_angle=d["flip_angle"], phase=d["phase"])
    else
        InstantGradient(; qvec=d["qvec"], origin=d["origin"], time=d["time"], apply_bvec=d["apply_bvec"])
    end
end

function parse_scanner(d::Dict)
    Scanner(
        B0=d["B0"],
        gradient=d["gradient"],
        slew_rate=d["slew_rate"],
    )
end

function parse_gradients(d::Dict)
    MRGradients(
        parse_shape(d["shape"], SVector{3, Float64}),
        d["origin"],
        d["apply_bvec"],
    )
end

function parse_shape(d::Dict, target_type=Float64)
    Shape(
        Float64.(d["times"]),
        target_type.(d["amplitudes"]),
    )
end

function parse_pulse(d::Dict)
    RFPulse(
        parse_shape(d["amplitude"]),
        parse_shape(d["phase"]),
        d["max_amplitude"]
    )
end

end