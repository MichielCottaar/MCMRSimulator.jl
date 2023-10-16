module Main
import StaticArrays: SVector
import Accessors: @set
import LinearAlgebra: ⋅
import ....Scanners: Scanner
import ...Methods: get_time, B0, get_rotation
import ..Gradients: MRGradients, gradient, rotate_bvec
import ..Instants: InstantComponent, Readout, InstantRFPulse, InstantGradient
import ..RadioFrequency: RFPulse, effective_pulse
import ..Shapes: ShapePart, sample_integral, sample
import ..Methods: start_time, end_time, add_TR
"""
    Sequence(;TR, components=[], scanner=Scanner(B0), B0=3.)

An MR sequence represented by a series of components repeated with a given repetition time (`TR`).

Possible sequence components are:
- [`RFPulse`](@ref): Radio-frequency pulse with user-provided amplitude and phase profile.
- [`InstantRFPulse`](@ref): instantaneous approximation of a radio-frequency pulse flipping the spin orientations.
- [`MRGradients`](@ref): MRI gradient profiles
- [`InstantGradient`](@ref): instantaneous gradients encoding spatial patterns in the spin phase distribution.
- [`Readout`](@ref): Store the spins at this timepoint.

The previous/current/next [`RFPulse`] at a specific time is given by [`previous_pulse`](@ref), [`current_pulse`](@ref), or [`next_pulse`](@ref). 
Equivalent functions are available to [`MRGradients`] called [`previous_gradient`](@ref), [`current_gradient`](@ref), or [`next_gradient`](@ref).
All of these will return `nothing` if there is no previous/current/next gradient/RF pulse.

The same functions exist for the previous/current/next instantaneous gradient/RF pulse (i.e., [`InstantRFPulse`](@ref) or [`InstantGradient`](@ref)), namely [`previous_instant`](@ref), [`current_instant`](@ref), or [`next_instant`](@ref). 

Note that all gradients/pulses repeat every `TR` milliseconds.
"""
struct Sequence{NI, NP, NR, NG}
    scanner :: Scanner
    gradients :: SVector{NG, MRGradients}
    instants :: SVector{NI, InstantComponent}
    pulses :: SVector{NP, RFPulse}
    TR :: Float64
    readout_times :: SVector{NR, Float64}
end

function Sequence(scanner::Scanner, components::AbstractVector, TR :: Real)
    @assert all([p isa Union{InstantComponent, Readout, RFPulse, MRGradients} for p in components])

    instants = sort([p for p in components if p isa InstantComponent], by=get_time)
    @assert length(instants) == 0 || (get_time(instants[end]) <= TR && get_time(instants[1]) >= 0)
    @assert length(instants) == length(unique(get_time.(instants)))

    readout_times = sort([get_time(p) for p in components if p isa Readout])
    @assert length(readout_times) == 0 || (readout_times[end] <= TR && readout_times[1] >= 0)

    gradients = sort([p for p in components if p isa MRGradients], by=start_time)
    for (p1, p2) in zip(gradients[1:end-1], gradients[2:end])
        @assert end_time(p1) <= start_time(p2)
    end
    for pulse in gradients
        @assert start_time(pulse) >= 0
        @assert end_time(pulse) <= TR
    end

    rf_pulses = sort([p for p in components if p isa RFPulse], by=start_time)
    for (p1, p2) in zip(rf_pulses[1:end-1], rf_pulses[2:end])
        @assert end_time(p1) <= start_time(p2)
    end
    for pulse in rf_pulses
        @assert start_time(pulse) >= 0
        @assert end_time(pulse) <= TR
    end
    Sequence{length(instants), length(rf_pulses), length(readout_times), length(gradients)}(scanner, gradients, instants, rf_pulses, Float64(TR), readout_times)
end

function Sequence(; scanner=nothing, components::AbstractVector=[], TR::Real, B0::Real=3.)
    if isnothing(scanner)
        scanner = Scanner(B0=B0)
    end
    Sequence(scanner, components, TR)
end

function Base.show(io::IO, seq::Sequence)
    if get(io, :compact, false)
        print(io, "Sequence($(length(seq.instants)) instants; $(length(seq.pulses)) RF pulses; $(length(seq.readouts)) readouts; TR=$(seq.TR)ms)")
    else
        text = "Sequence (TR=$(seq.TR)ms):\n"
        for pulse in sort([seq.pulses..., seq.instants..., seq.readout_times...], by=start_time)
            if pulse isa Float64
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

start_time(pulse::Union{InstantComponent, Readout}) = pulse.time
end_time(pulse::Union{InstantComponent, Readout}) = pulse.time
start_time(number::Number) = number
end_time(number::Number) = number

function add_TR(pulse::Union{InstantComponent, Readout}, delta_time) 
    @set pulse.time += delta_time
end

next_index(pulses, time) = searchsortedfirst(pulses, time, by=end_time)
previous_index(pulses, time) = searchsortedlast(pulses, time, by=start_time)
function current_index(pulses, time)
    index = next_index(pulses, time)
    return index == previous_index(pulses, time) ? index : -1
end

for relative in ("previous", "current", "next")
    for (pulse_type, type_name) in [
        ("pulse", "[`RFPulse`](@ref)")
        ("gradient", "[`MRGradients`](@ref)")
        ("instant", "[`InstantGradient`](@ref) or [`InstantRFPulse`](@ref)")
    ]
        func_name = Symbol(relative * "_" * pulse_type)
        plural_pulse_type = Symbol(pulse_type * "s")
        index_func = Symbol(relative * "_index")
        @eval begin
            """
                $($func_name)(sequence, time)
            
            Returns the $($relative) $($type_name) in the [`Sequence`](@ref) relative to `time`.
            If there is no $($relative) $($pulse_type) `nothing` is returned instead.
            """
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
    gradient([position, ], sequence, t1)

Gets the off-resonance field (units: kHz) at given `position` at time `t1`.
If position is not supplied the MR gradient is returned as a vector (units: kHz/um).
"""
function gradient(position::AbstractVector, seq::Sequence, time::Number) 
    grad = current_gradient(seq, time)
    isnothing(grad) ? zero(Float64) : gradient(position, current_gradient(seq, time), mod(time, seq.TR))
end

function gradient(seq::Sequence, time::Number)
    grad = current_gradient(seq, time)
    isnothing(grad) ? zero(Float64) : gradient(current_gradient(seq, time), mod(time, seq.TR))
end

"""
    SequencePart(sequence, t1, t2)

Represents a small part of a [`Sequence`](@ref) between times `t1` and `t2`.
During this time the RF pulse and gradients are assumed to change linearly and hence can be represented as [`ShapePart`](@ref) objects.
This is used to guide the spin evolution during the timestep between times `t1` and `t2` in the simulation.
"""
struct SequencePart
    rf_amplitude::ShapePart{Float64}
    rf_phase::ShapePart{Float64}
    gradients::ShapePart{SVector{3, Float64}}
    origin::SVector{3, Float64}
    total_time::Float64
    B0::Float64
end

function SequencePart(sequence::Sequence, t1::Number, t2::Number)
    nTR = Int(div((t1 + t2) / 2, sequence.TR, RoundDown))
    t1_norm = t1 - nTR * sequence.TR
    t2_norm = t2 - nTR * sequence.TR
    pulse = current_pulse(sequence, (t1_norm + t2_norm) / 2)
    if isnothing(pulse)
        amp = phase = ShapePart(zero(Float64), zero(Float64), zero(Float64))
    else
        amp = ShapePart(pulse.amplitude, t1_norm, t2_norm)
        phase = ShapePart(pulse.phase, t1_norm, t2_norm)
    end
    grad = current_gradient(sequence, (t1_norm + t2_norm) / 2)
    if isnothing(grad)
        grad_part = ShapePart(zero(SVector{3, Float64}), zero(SVector{3, Float64}), zero(SVector{3, Float64}))
    else
        grad_part = ShapePart(grad.shape, t1_norm, t2_norm)
    end
    SequencePart(
        amp, phase,
        grad_part, isnothing(grad) ? zero(SVector{3, Float64}) : grad.origin,
        t2 - t1, B0(sequence)
    )
end

gradient(position, part::SequencePart, t0, t1) = sample(part.gradients, t0, t1) ⋅ (position - part.origin)

effective_pulse(part::SequencePart, t0::Number, t1::Number) = InstantRFPulse(flip_angle=sample_integral(part.rf_amplitude, t0, t1) * part.total_time * 360., phase=sample(part.rf_phase, t0, t1))
B0(part::SequencePart) = part.B0

"""
    rotate_bvec(sequence, bvec)

Rotates diffusion-weighting gradients in `sequence` by `bvec`.
The rotation is only applied to [`MRGradients`](@ref) and [`InstantGradient`](@ref) objects in the sequence for which `apply_bvec` is true.
This should include any diffusion-weighted gradients, but not slice-selection, crusher, or readout gradients.

`bvec` can be a full rotation matrix or a vector.
In the latter case the x-axis of the gradients will be rotated to the vector using a minimal rotation (see [`get_rotation`](@ref)).
"""
function rotate_bvec(sequence::Sequence{I, P, R, G}, bvec) where {I, P, R, G}
    rotation = get_rotation(bvec, 3)
    function apply_rotation(grad::MRGradients)
        return grad.apply_bvec ? rotate_bvec(grad, rotation) : grad
    end
    function apply_rotation(instant_grad::InstantGradient)
        return instant_grad.apply_bvec ? rotate_bvec(instant_grad, rotation) : instant_grad
    end
    function apply_rotation(pulse::InstantRFPulse)
        return pulse
    end
    return Sequence{I, P, R, G}(
        sequence.scanner,
        apply_rotation.(sequence.gradients),
        apply_rotation.(sequence.instants),
        sequence.pulses,
        sequence.TR,
        sequence.readout_times
    )
end


"""
    can_rotate_bvec(sequence)

Checks whether [`rotate_bvec`](@ref) will have any effect on the sequence (i.e., `apply_bvec` has been set to true for any non-zero gradients in the sequence).
"""
function can_rotate_bvec(sequence::Sequence)
    can_rotate(grad::MRGradients) = grad.apply_bvec && any(norm.(grad.shape.amplitude) .> 0)
    can_rotate(grad::InstantGradient) = grad.apply_bvec && norm(grad.qvec) > 0
    can_rotate(no_grad::InstantRFPulse) = false
    return any(can_rotate.(sequence.gradients)) || any(can_rotate.(sequence.instants))
end

end