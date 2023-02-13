include("instants.jl")
include("shape.jl")
include("gradients.jl")
include("radio_frequency.jl")

"""
    Sequence(;TR, gradients=nothing, pulses=nothing, scanner=Scanner(B0), B0=3., interplate_gradients=:step)

An MR sequence represented by a series of pulses repeated with a given repetition time (`TR`).

Possible sequence components are:
- [`RFPulse`](@ref): Radio-frequency pulse with user-provided amplitude and phase profile.
- [`InstantRFPulse`](@ref): instantaneous approximation of a radio-frequency pulse flipping the spin orientations.
- [`InstantGradient`](@ref): instantaneous gradients encoding spatial patterns in the spin phase distribution.
- [`Readout`](@ref): Store the spins at this timepoint.

The previous/current/next RF pulse at a specific time is given by [`previous_pulse`](@ref), [`current_pulse`](@ref), or [`next_pulse`](@ref). 
All of these will return `nothing` if there is no previous/current/next pulse.
The same functions exist for the previous/current/next instantaneous pulse (i.e., [`InstantRFPulse`](@ref) or [`InstantGradient`](@ref)), 
    namely [`previous_instant`](@ref), [`current_instant`](@ref), or [`next_instant`](@ref). 
Note that all gradients/pulses repeat every `TR`.
"""
struct Sequence{NI, NP, NR, G<:MRGradients}
    scanner :: Scanner
    gradient :: G
    instants :: SVector{NI, InstantComponent}
    pulses :: SVector{NP, RFPulse}
    TR :: Float
    readout_times :: SVector{NR, Float}
    function Sequence(scanner::Scanner, gradients::MRGradients, pulses::AbstractVector, TR :: Real)
        all([p isa Union{InstantComponent, Readout, RFPulse} for p in pulses])
        instants = sort([p for p in pulses if p isa InstantComponent], by=get_time)
        @assert length(instants) == 0 || (get_time(instants[end]) <= TR && get_time(instants[1]) >= 0)
        @assert length(instants) == length(unique(get_time.(instants)))
        readout_times = sort([get_time(p) for p in pulses if p isa Readout])
        @assert length(readout_times) == 0 || (readout_times[end] <= TR && readout_times[1] >= 0)
        rf_pulses = sort([p for p in pulses if p isa RFPulse], by=start_time)
        for (p1, p2) in zip(rf_pulses[1:end-1], rf_pulses[2:end])
            @assert end_time(p1) <= start_time(p2)
        end
        for pulse in rf_pulses
            @assert start_time(pulse) >= 0
            @assert end_time(pulse) <= TR
        end
        new{length(instants), length(rf_pulses), length(readout_times), typeof(gradients)}(scanner, gradients, instants, rf_pulses, Float(TR), readout_times)
    end
end

function Sequence(; scanner=nothing, gradients=nothing, pulses::AbstractVector=[], TR::Real, B0::Real=3.)
    if isnothing(scanner)
        scanner = Scanner(B0=B0)
    end
    if isnothing(gradients)
        gradients = MRGradients()
    elseif !isa(gradients, MRGradients)
        gradients = MRGradients(gradients)
    end
    Sequence(scanner, gradients, pulses, TR)
end

function Base.show(io::IO, seq::Sequence)
    if get(io, :compact, false)
        print(io, "Sequence($(length(seq.instants)) instants; $(length(seq.pulses)) RF pulses; $(length(seq.readouts)) readouts; TR=$(seq.TR)ms)")
    else
        text = "Sequence (TR=$(seq.TR)ms):\n"
        for pulse in sort([seq.pulses..., seq.instants..., seq.readout_times...], by=start_time)
            if pulse isa Float
                text *= "    - Readout at $(pulse)ms\n"
            else
                text *= "    - $(pulse)\n"
            end
        end
        print(io, text)
    end
end

"""
    effective_pulse(sequence, t1, t2)

Returns the [`InstantRFPulse`](@ref) that has the same effect as the radio-frequency pulses ([`RFPulse`](@ref)) in the provided sequence will have between `t1` and `t2`.
"""
function effective_pulse(sequence::Sequence, t1::Number, t2::Number)
    current = current_pulse(sequence, (t1 + t2) / 2)
    if isnothing(current)
        return InstantRFPulse(0, 0, 0)
    else
        @assert t1 >= start_time(current)
        @assert t2 <= end_time(current)
        return effective_pulse(current, t1, t2)
    end
end

Scanner(sequence::Sequence) = sequence.scanner
B0(sequence::Sequence) = B0(Scanner(sequence))

start_time(pulse::InstantComponent) = pulse.time
end_time(pulse::InstantComponent) = pulse.time
start_time(number::Number) = number
end_time(number::Number) = number

function add_TR(pulse::InstantComponent, delta_time) 
    @set pulse.time += delta_time
end

next_index(pulses, time) = searchsortedfirst(pulses, time, by=end_time)
previous_index(pulses, time) = searchsortedlast(pulses, time, by=start_time)
function current_index(pulses, time)
    index = next_index(pulses, time)
    return index == previous_index(pulses, time) ? index : -1
end

for relative in ("previous", "current", "next")
    for pulse_type in ("pulse", "instant")
        func_name = Symbol(relative * "_" * pulse_type)
        plural_pulse_type = Symbol(pulse_type * "s")
        index_func = Symbol(relative * "_index")
        @eval begin
            function $func_name(sequence, time)
                nTR = Int(div(time, sequence.TR, RoundDown))
                within_TR = time - sequence.TR * nTR
                all_pulses = sequence.$plural_pulse_type
                if length(all_pulses) == 0
                    return nothing
                end
                index = $index_func(all_pulses, within_TR)
                if index == 0
                    nTR -= 1
                    index = length(all_pulses)
                elseif index == length(all_pulses) + 1
                    nTR += 1
                    index = 1
                end
                if index == -1 || nTR < 0
                    return nothing
                end
                add_TR(all_pulses[index], nTR * sequence.TR)
            end
        end
    end
end

"""
    gradient([position, ], sequence, t1[, t2])

Gets the off-resonance field (units: kHz) at given `position` at time `t1` (or integrated between times `t1` and `t2`).
If position is not supplied the MR gradient is returned as a [`PosVector`](@ref) (units: kHz/um).
"""
gradient(position::AbstractVector, seq::Sequence, time::Number) = gradient(position, seq.gradient, mod(time, seq.TR))
gradient(position::AbstractVector, seq::Sequence, time1::Number, time2::Number) = gradient(position, seq.gradient, mod(time1, seq.TR), iszero(mod(time2, seq.TR)) ? seq.TR : mod(time2, seq.TR))
gradient(seq::Sequence, time::Number) = gradient(seq.gradient, mod(time, seq.TR))
gradient(seq::Sequence, time1::Number, time2::Number) = gradient(seq.gradient, mod(time1, seq.TR), iszero(mod(time2, seq.TR)) ? seq.TR : mod(time2, seq.TR))

"""
    SequencePart(sequence, t1, t2)

Represents a small part of a [`Sequence`](@ref) between times `t1` and `t2`.
During this time the RF pulse and gradients are assumed to change linearly and hence can be represented as [`ShapePart`](@ref) objects.
This is used to guide the spin evolution during the timestep between times `t1` and `t2` during the simulation.
"""
struct SequencePart
    rf_amplitude::ShapePart
    rf_phase::ShapePart
    Gx::ShapePart
    Gy::ShapePart
    Gz::ShapePart
    origin::PosVector
    total_time::Float
    B0::Float
end

function SequencePart(sequence::Sequence, t1::Number, t2::Number)
    t1_norm = mod(t1, sequence.TR)
    t2_norm = mod(t2, sequence.TR)
    pulse = current_pulse(sequence, (t1_norm + t2_norm) / 2)
    if isnothing(pulse)
        amp = phase = ShapePart(zero(Float), zero(Float), zero(Float))
    else
        amp = ShapePart(pulse.amplitude, t1_norm, t2_norm)
        phase = ShapePart(pulse.phase, t1_norm, t2_norm)
    end
    SequencePart(
        amp, phase,
        ShapePart(sequence.gradient.Gx, t1_norm, t2_norm),
        ShapePart(sequence.gradient.Gy, t1_norm, t2_norm),
        ShapePart(sequence.gradient.Gz, t1_norm, t2_norm),
        sequence.gradient.origin,
        t2 - t1,
        B0(sequence)
    )
end

function gradient(position, part::SequencePart, t1, t2)
    grad = zero(Float)
    for (index, shape_part) in zip(1:3, [part.Gx, part.Gy, part.Gz])
        grad += amplitude(shape_part, t1, t2) * (position[index] - part.origin[index])
    end
    grad
end

effective_pulse(part::SequencePart, t0::Number, t1::Number) = InstantRFPulse(flip_angle=amplitude(part.rf_amplitude, t0, t1) * (t1 - t0) * part.total_time * 360., phase=amplitude(part.rf_phase, t0, t1))
B0(part::SequencePart) = part.B0

include("diffusion.jl")
include("pulseseq.jl")