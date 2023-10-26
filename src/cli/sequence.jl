"""
Defines the command line interface to `mcmr sequence`.
"""
module Sequence

import ArgParse: ArgParseSettings, @add_arg_table!, add_arg_table!, add_arg_group!, parse_args
import ...SequenceBuilder.Sequences.SpinEcho: spin_echo, dwi
import ...SequenceBuilder.Sequences.GradientEcho: gradient_echo
import ...SequenceBuilder.Diffusion: gen_crusher
import ...Sequences: InstantRFPulse, constant_pulse, write_sequence
import ...Scanners: Scanner, predefined_scanners, max_gradient


function known_sequence_parser(name; kwargs...)
    parser = ArgParseSettings(prog="mcmr sequence $name"; kwargs...)
    @add_arg_table! parser begin
        "output-file"
            help = "Create a sequence JSON file containing the $name sequence"
            required = true
    end

    add_arg_group!(parser, "Scanner parameters", :scanner)
    @add_arg_table! parser begin
        "--scanner"
            help = "predefined scanners. One of $(join(keys(predefined_scanners), ", ", ", or ")). If set all other scanner parameter are ignored."
            range_tester = (x -> x in keys(predefined_scanners))
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
            default = Inf
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
    scanner_name = pop!(arguments, "scanner")
    if ~isnothing(scanner_name)
        for (field, ref_value) in [("kHz", false), ("B0", 3.), ("max-gradient", Inf), ("max-slew-rate", Inf)]
            value = pop!(arguments, field)
            if value != ref_value
                @warn "--$field flag will be ignored, because --scanner is set."
            end
        end
        return predefined_scanners[scanner_name]
    end
    units = pop!(arguments, "kHz") ? :kHz : :Tesla
    B0 = pop!(arguments, "B0")
    grad = pop!(arguments, "max-gradient")
    slew = pop!(arguments, "max-slew-rate")
    return Scanner(B0=B0, gradient=grad, slew_rate=slew, units=units)
end

function add_pulse_to_parser!(parser, pulse_name; flip_angle=90., phase=0., duration=0., with_crusher=false, crusher_qval=nothing, crusher_duration=nothing)
    if iszero(length(pulse_name))
        addition = "-"
    else
        addition = "--$pulse_name"
    end
    caps = uppercasefirst(pulse_name)
    add_arg_table!(parser,
        "$(addition)-duration", Dict(
            :help => "$caps pulse duration (ms).",
            :arg_type => Float64,
            :default => Float64(duration),
        ),
        "$(addition)-flip-angle", Dict(
            :help => "$caps pulse flip angle (deg).",
            :arg_type => Float64,
            :default => Float64(flip_angle),
        ),
        "$(addition)-phase", Dict(
            :help => "$caps pulse phase (deg).",
            :arg_type => Float64,
            :default => Float64(phase),
        ),
    )
    if with_crusher
        add_crusher_to_parser!(parser, pulse_name * "-crusher"; qval=crusher_qval, duration=crusher_duration, description="crusher around $caps pulse")
    end
end

function add_crusher_to_parser!(parser, crusher_name; qval=nothing, duration=nothing, description=nothing)
    if isnothing(description)
        description = crusher_name * " gradients"
    end
    if iszero(length(crusher_name))
        addition = "-"
    else
        addition = "--$crusher_name"
    end
    if isnothing(qval)
        qval = Inf
    end
    if isnothing(duration)
        duration = Inf
    end
    add_arg_table!(parser,
        "$(addition)-qval", Dict(
            :help => "Strength of $description (rad/um).",
            :arg_type => Float64,
            :default => Float64(qval),
        ),
        "$(addition)-duration", Dict(
            :help => "Duration of $description (ms).",
            :arg_type => Float64,
            :default => Float64(duration),
        ),
    )
end

function get_crusher(arguments, crusher_name, scanner; default_qval=1.)
    addition = iszero(length(crusher_name)) ? "" : "$(crusher_name)-"

    qval = pop!(arguments, "$(addition)qval")
    duration = pop!(arguments, "$(addition)duration")
    @show crusher_name qval duration
    if isinf(qval) && (isinf(duration) || isinf(max_gradient(scanner)))
        qval = default_qval
        if isinf(duration)
            duration = 0.
        end
    end
    if iszero(qval)
        println(duration)
        return isinf(duration) ? 0. : duration
    end

    if isinf(qval)
        qval = nothing
    end
    if isinf(duration)
        duration = nothing
    end
    result = gen_crusher(qval=qval, duration=duration, scanner=scanner)
    @show result
    return result
end

function get_pulse(arguments, pulse_name, scanner::Scanner)
    addition = iszero(length(pulse_name)) ? "" : "$(pulse_name)-"
    duration = pop!(arguments, "$(addition)duration")
    fa = pop!(arguments, "$(addition)flip-angle")
    phase = pop!(arguments, "$(addition)phase")
    if iszero(duration)
        pulse = InstantRFPulse(flip_angle=fa, phase=phase)
    else
        pulse = constant_pulse(duration, fa; phase0=phase)
    end
    if "$(addition)crusher-qval" in keys(arguments)
        crusher = get_crusher(arguments, "$(pulse_name)-crusher", scanner; default_qval=0.)
        return [crusher, pulse, crusher]
    else
        return pulse
    end
end


function run_dwi(args=ARGS::AbstractVector[<:AbstractString]; kwargs...)
    parser = known_sequence_parser("dw-pgse"; kwargs...)
    parser.description = "Implement diffusion-weighted (DW) pulsed-gradient spin-echo (PGSE) sequence."

    @add_arg_table! parser begin
        "--TE"
            help = "Spin echo and readout time in ms."
            arg_type = Float64
            required = true
        "-b", "--bval"
            help = "Strength of the diffusion-weighting (ms/um^2)."
            arg_type = Float64
            required = true
        "--gradient-duration"
            help = "Duration of the gradients (ms). Default: maximum possible"
            arg_type = Float64
        "--readout-time"
            help = "Duration of the readout (ms). The actual signal readout will happen half-way this period."
            arg_type = Float64
            default = 0.
    end
    add_pulse_to_parser!(parser, "excitation"; flip_angle=90, phase=-90, duration=0)
    add_pulse_to_parser!(parser, "refocus"; flip_angle=180, phase=0, duration=0, with_crusher=true)
    add_crusher_to_parser!(parser, "crusher"; description="crusher gradient after readout")


    as_dict = parse_args(args, parser)
    output_file = pop!(as_dict, "output-file")
    scanner = get_scanner(as_dict)
    as_dict["scanner"] = scanner
    as_dict["excitation_pulse"] = get_pulse(as_dict, "excitation", scanner)
    as_dict["refocus_pulse"] = get_pulse(as_dict, "refocus", scanner)
    as_dict["crusher"] = get_crusher(as_dict, "crusher", scanner)
    sequence = dwi(; Dict(Symbol(replace(k, "-"=>"_")) => v for (k, v) in as_dict)...)
    write_sequence(output_file, sequence)
end


function run_spin_echo(args=ARGS::AbstractVector[<:AbstractString]; kwargs...)
    parser = known_sequence_parser("spin_echo"; kwargs...)
    parser.description = "Implement spin echo sequence with single readout."

    @add_arg_table! parser begin
        "--TE"
            help = "Spin echo time in ms."
            arg_type = Float64
            required = true
        "--readout-time"
            help = "Duration of the readout (ms). The actual signal readout will happen half-way this period."
            arg_type = Float64
            default = 0.
    end
    add_pulse_to_parser!(parser, "excitation"; flip_angle=90, phase=-90, duration=0)
    add_pulse_to_parser!(parser, "refocus"; flip_angle=180, phase=0, duration=0)
    add_crusher_to_parser!(parser, "crusher"; description="crusher gradient after readout")


    as_dict = parse_args(args, parser)
    output_file = pop!(as_dict, "output-file")
    scanner = get_scanner(as_dict)
    as_dict["scanner"] = scanner
    as_dict["excitation_pulse"] = get_pulse(as_dict, "excitation", scanner)
    as_dict["refocus_pulse"] = get_pulse(as_dict, "refocus", scanner)
    as_dict["crusher"] = get_crusher(as_dict, "crusher", scanner)
    TE = pop!(as_dict, "TE")
    sequence = spin_echo(TE; Dict(Symbol(replace(k, "-"=>"_")) => v for (k, v) in as_dict)...)
    write_sequence(output_file, sequence)
end


function run_gradient_echo(args=ARGS::AbstractVector[<:AbstractString]; kwargs...)
    parser = known_sequence_parser("gradient_echo"; kwargs...)
    parser.description = "Implement gradient echo sequence with single readout."

    @add_arg_table! parser begin
        "--TE"
            help = "Gradient echo time in ms."
            arg_type = Float64
            required = true
    end
    add_pulse_to_parser!(parser, "excitation"; flip_angle=90, phase=-90, duration=0)
    add_crusher_to_parser!(parser, "crusher"; description="crusher gradient after readout")


    as_dict = parse_args(args, parser)
    output_file = pop!(as_dict, "output-file")
    scanner = get_scanner(as_dict)
    as_dict["scanner"] = scanner
    as_dict["excitation_pulse"] = get_pulse(as_dict, "excitation", scanner)
    as_dict["crusher"] = get_crusher(as_dict, "crusher", scanner)
    TE = pop!(as_dict, "TE")
    sequence = gradient_echo(TE; Dict(Symbol(replace(k, "-"=>"_")) => v for (k, v) in as_dict)...)
    write_sequence(output_file, sequence)
end

"""
    run_main([arguments])

Runs the `mcmr sequence` command line interface.
Arguments are provided as a sequence of strings.
By default it is set to `ARGS`.
"""
function run_main(args=ARGS::AbstractVector[<:AbstractString]; kwargs...)
    pre_created = Dict(
        "dwi" => run_dwi,
        "dw-pgse" => run_dwi,
        "spin-echo" => run_spin_echo,
        "gradient-echo" => run_gradient_echo,
    )
    if length(args) == 0
        println(stderr, "No mcmr sequence sub-command given.\n")
    else
        if haskey(pre_created, args[1])
            return pre_created[args[1]](args[2:end]; kwargs...)
        else
            println(stderr, "Invalid mcmr command sequence  $(args[1]) given.\n")
        end
    end
    names = join(keys(pre_created), "/")
    println(stderr, "usage: mcmr sequence {$names}")
    return Cint(1)

end

end