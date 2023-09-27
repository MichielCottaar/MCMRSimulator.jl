"""
Defines the command line interface to `mcmr sequence`.
"""
module Sequence
import ArgParse: ArgParseSettings, @add_arg_table!, add_arg_group!, parse_args
import ...SequenceBuilder.Sequences.SpinEcho: spin_echo, dwi
import ...Sequences: InstantRFPulse, constant_pulse
import ...Scanners: Scanner


function known_sequence_parser(name)
    parser = ArgParseSettings(prog="mcmr sequence $name")
    @add_arg_table! parser begin
        "output-file"
            help = "Create a sequence JSON file containing the $name sequence"
    end

    add_arg_group!(parser, "Scanner parameters", :scanner)
    @add_arg_table! parser begin
        "--B0"
            help = "Magnetic field strength in Tesla."
            arg_type = Float64
            default = 3.
        "--max-gradient"
            help = "Gradient strength in mT/m (unless --kHz is set)."
            arg_type = Float64
            default = Inf
        "--max-slew-rate"
            help = "Slew rate in T/m/s (unless --kHz is set)."
            arg_type = Float64
            default = 0.
        "--kHz"
            help = "Give gradient in kHz/um and slew rate in kHz/um/ms rather than the default values"
            action = :store_true
    end

    add_arg_group!(parser, "Sequence-specific parameters", :sequence)
    @add_arg_table! parser begin
        "--TR"
            help = "Sequence repetition time (ms)."
            arg_type = Float64
            required = true
    end
end

function get_scanner(arguments)
    units = pop!(parser, "kHz") ? :kHz ? :Tesla
    B0 = pop!(parser, "B0")
    grad = pop!(parser, "max-gradient")
    slew = pop!(parser, "max-slew-rate")
    return Scanner(B0=B0, gradient=grad, slew_rate=slew, units=units)
end

function add_pulse_to_parser!(parser, pulse_name; flip_angle=90., phase=0., duration=0.)
    if iszero(length(pulse_name))
        addition = "-"
    else
        addition = "--$pulse_name"
    end
    caps = uppercasefirst(pulse_name)
    @add_arg_table! parser begin
        "$(addition)-duration"
            help = "$caps pulse duration (ms)."
            arg_type = Float64
            default = duration
        "$(addition)-flip-angle"
            help = "$caps pulse flip angle (deg)."
            arg_type = Float64
            default = flip_angle
        "$(addition)-phase"
            help = "$caps phase (deg)."
            arg_type = Float64
            default = phase
    end
end

function get_pulse(arguments, pulse_name)
    addition = iszero(length(pulse_name)) ? "" : "$(pulse_name)-"
    duration = pop!(arguments, "$(addition)duration")
    fa = pop!(arguments, "$(addition)flip-angle")
    phase = pop!(arguments, "$(addition)phase")
    if iszero(duration)
        return InstantRFPulse(flip_angle=fa, phase=phase)
    else
        return constant_pulse(duration, fa; phase0=phase)
    end
end


function run_dwi(args=ARGS::AbstractVector[<:AbstractString])
    parser = known_sequence_parser("dwi")

    @add_arg_table! parser begin
        "--TE"
            help = "Echo time in ms."
            arg_type = Float64
            required = true
        "--bval"
            help = "Strength of the diffusion-weighting (ms/um^2)."
            arg_type = Float64
        "--gradient-duration"
            help = "Duration of the gradients (ms). Default: maximum possible"
            arg_type = Float64
        "--readout-time"
            help = "Duration of the readout (ms). The actual signal readout will happen half-way this period."
            arg_type = Float64
            default = 0.
    end
    add_pulse_to_parser!(parser, "excitation"; flip_angle=90, phase=-90, duration=0)
    add_pulse_to_parser!(parser, "refocus"; flip_angle=180, phase=0, duration=0)


    as_dict = parse_args(args)
    output_file = pop!(as_dict, "output-file")
    as_dict["excitation_pulse"] = get_pulse(as_dict, "excitation")
    as_dict["refocus_pulse"] = get_pulse(as_dict, "refocus")
    as_dict["scanner"] = get_scanner(as_dict)
    sequence = dwi(; as_dict...)
    write_sequence()
end


"""
    run_main([arguments])

Runs the `mcmr sequence` command line interface.
Arguments are provided as a sequence of strings.
By default it is set to `ARGS`.
"""
function run_main(args=ARGS::AbstractVector[<:AbstractString])
    if length(args) == 0
        println("No mcmr sequence sub-command given.\n")
    else
        if args[1] == "custom"
            return run_custom(args[2:end])
        elseif args[1] == "dwi"
            return Geometry.run_dwi(args[2:end])
        else
            println("Invalid mcmr command sequence  $(args[1]) given.\n")
        end
    end
    println("usage: mcmr {custom/}")
    return Cint(1)

end

end