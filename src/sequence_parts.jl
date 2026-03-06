"""
Breaks up one or more `MRIBuilder.Sequence` objects into individual timesteps.

The main interface is by calling `iter_parts`, which will result into a vector of some subtypes of `SequencePart`.

To just get the readouts call [`MCMRSimulator.get_readouts`](@ref MCMRSimulator.SequenceParts.get_readouts).
"""
module SequenceParts
import StaticArrays: SVector
import LinearAlgebra: norm
import Statistics: mean
import Interpolations: linear_interpolation, deduplicate_knots!
import ..TimeSteps: TimeStep
import ..Spins: static_vector_type
import KomaMRIBase
import KomaMRIFiles: read_seq


"""
A short part of the sequence that can be handled by the simulator.

There are 4 types:
- [`EmptyPart`](@ref): no pulse or gradient
- [`ContantPart`](@ref): no pulse; constant gradient
- [`LinearPart`](@ref): no pulse; linear gradient
- [`PulsePart`](@ref): RF pulse; gradient can be any of the above
"""
abstract type SequencePart end

"""
Any part of the sequence not containing an RF pulse.
"""
abstract type NoPulsePart <: SequencePart end

"""
Part of the sequence with no RF pulse or gradient.
"""
struct EmptyPart <: NoPulsePart
end


"""
Part of the sequence with no RF pulse and a constant gradient amplitude (given by `strength`).
"""
struct ConstantPart <: NoPulsePart
    strength :: SVector{3, Float64}
end


"""
Part of the sequence with no RF pulse and a linearly changing gradient amplitude (given by `start` and `final`).
"""
struct LinearPart <: NoPulsePart
    start :: SVector{3, Float64}
    final :: SVector{3, Float64}
end


"""
    ConstantPulse(amplitude, phase, frequency)

Stores the RF pulse state during a single timestep.

It contains:
- `amplitude` of RF pulse in kHz (assumed to be constant over the timestep).
- `phase` of the RF pulse at the beginning of the timestep in degrees.
- off-resonance `frequency` of the RF pulse in kHz.
"""
struct ConstantPulse
    amplitude :: Float64
    phase :: Float64
    frequency :: Float64
end

"""
    PulsePart{T}(pulse, gradient)

A `SequencePart` containing a constant RF pulse and a gradient of type `T` (which is a subtype of `NoPulsePart`).

The RF pulse is represented by a vector of `ConstantPulse`s, which are assumed to be played sequentially during the part.
The duration of the first and last `ConstantPulse` might be different from the others as represented by `first_duration` and `last_duration`.
If `first_duration` and `last_duration` are 1 (default) then all `ConstantPulse`s are assumed to have the same duration.

In the special case of there only being one `ConstantPulse`, only `first_duration` is used and represents the duration of that `ConstantPulse`.
"""
struct PulsePart{T<:NoPulsePart} <: SequencePart
    pulse :: Vector{ConstantPulse}
    gradient :: T
    first_duration :: Float64
    last_duration :: Float64
end

PulsePart(pulse::Vector{ConstantPulse}, gradient::T) where {T<:NoPulsePart} = PulsePart{T}(pulse, gradient, 1., 1.)


function split_fraction(part::PulsePart{T}, fraction::Number) where {T<:NoPulsePart}
    grad1, grad2 = split_fraction(part.gradient, fraction)
    total_duration = length(part.pulse) - 2 + part.first_duration + part.last_duration
    split_duration = total_duration * fraction

    if split_duration < part.first_duration
        # split is within the first pulse
        return (
            PulsePart(part.pulse[1:1], grad1, split_duration, 1.),
            PulsePart(part.pulse, grad2, part.first_duration - split_duration, part.last_duration)
        )
    elseif split_duration > total_duration - part.last_duration
        # split is within the last pulse
        return (
            PulsePart(part.pulse, grad1, part.first_duration, part.last_duration - (total_duration - split_duration)),
            PulsePart(part.pulse[end:end], grad2, 1., split_duration - (total_duration - part.last_duration))
        )
    end
    index_split = split_duration - part.first_duration + 1
    if index_split == round(index_split)
        # split is between two pulses
        cut_index = Int(index_split)
        return (
            PulsePart(part.pulse[1:cut_index], grad1, part.first_duration, 1.),
            PulsePart(part.pulse[cut_index+1:end], grad2, 1., part.last_duration)
        )
    end

    cut_index = Int(ceil(index_split))
    relative_split = index_split - (cut_index - 1)
    return (
        PulsePart(part.pulse[1:cut_index], grad1, part.first_duration, relative_split),
        PulsePart(part.pulse[cut_index:end], grad2, 1 - relative_split, part.last_duration)
    )
end


"""
Instantaneous event within the sequence, such as a readout or an narrow-pulse approximation of a gradient or RF pulse.
"""
abstract type SequenceEvent end

"""
Instantaneous approximation of an RF pulse with the `flip_angle` in degrees and `phase` in degrees.

The bandwidth is infinite, so the pulse does not have a frequency.
"""
struct PulseEvent <: SequenceEvent
    flip_angle :: Float64
    phase :: Float64
end

"""
Instantaneous approximation of a gradient with the q-vector `qvec` in 1/um.
"""
struct GradientEvent <: SequenceEvent
    qvec :: SVector{3, Float64}
end

"""
    IndexedReadout(time, TR, readout)

Represents a readout in an MR sequence.

Generate these for a sequence of interest using [`get_readouts`](@ref).

All indices are integers. They refer to:
- `time`: time since beginning of this repetition of the sequence in ms.
- `TR`: which TR the simulation is in (defaults to 0 if not set).
- `readout`: which readout within the TR (or total sequence) this is.
"""
struct IndexedReadout <: SequenceEvent
    time::Float64
    TR::Int
    readout::Int
end


"""
    SequenceWaveform(gradients, rf_pulses, samples, TR)
    SequenceWaveform(koma_sequence::KomaMRIBase.Sequence)
    SequenceWaveform(pulseq_filename::String)

Represents the sequence as a continous gradient waveform, a list of RF pulses and the readout times.
"""
struct SequenceWaveform
    grads::NTuple{3, Tuple{Vector{Float64}, Vector{Float64}}}
    rf::Vector{Tuple{Float64, Float64, Vector{ConstantPulse}}}
    samples::Vector{Float64}
    TR::Float64
end

function SequenceWaveform(sequence::KomaMRIBase.Sequence)
    grads = (
        gradient_waveform(sequence, 1),
        gradient_waveform(sequence, 2),
        gradient_waveform(sequence, 3),
    )
    rf = get_pulses(sequence)
    samples = readout_times(sequence)
    return SequenceWaveform(
        grads, 
        rf, 
        samples,
        grads[1][1][end]
    )
end

SequenceWaveform(filename::AbstractString) = SequenceWaveform(read_seq(filename))


"""
    parts(sequences::AbstractVector, start_time, timestep; kwargs...)

Splits the sequences into parts that can be handled by the simulator based on the control times of the sequence and the maximum allowed timestep.

The control times include:
- any readout times. The readout times are determined by `get_readouts` and can be altered by passing on the relevant keywords to `parts`.
- any time points where the gradient slew rate changes.
- any time points where an RF pulse starts or ends.
"""
function parts(sequences::AbstractVector, start_time::Number, timestep::TimeStep; kwargs...)
    waveforms = [SequenceWaveform(seq) for seq in sequences]
    (repeats, readouts) = zip([get_readouts(waveform.samples, waveform.TR, start_time; kwargs...) for waveform in waveforms]...)

    # Figure out control times
    set_control_times = Set{Float64}([start_time])
    for ro in readouts
        union!(set_control_times, [r.time for r in ro])
    end
    max_ro = maximum(set_control_times)
    for waveform in waveforms
        for index in 1:3
            (new_times, _) = compress_timeseries(waveform.grads[index]...)

            union!(set_control_times, filter(new_times) do t
                t > start_time && t < max_ro
            end)
        end
    end

    for (t0_rf, t1_rf, _) in waveforms[1].rf
        for t in (t0_rf, t1_rf)
            if start_time < t < max_ro
                push!(set_control_times, t)
            end
        end
    end
    control_times = sort(collect(set_control_times))

    grad_interpolators = map(waveforms, repeats) do waveform, repeating
        map(waveform.grads) do (times, ampls)
            inter = linear_interpolation(deduplicate_knots!(times), ampls; extrapolation_bc=0.)
            if repeating
                time -> inter(time % waveform.TR)
            else
                inter
            end
        end
    end

    # split control times based on the maximum allowed timestep
    for (t0, t1) in zip(control_times[1:end-1], control_times[2:end])
        max_gradient = maximum(map(grad_interpolators) do interpolators
            maximum(map(interpolators) do interpolator
                max(abs(interpolator(t0)), abs(interpolator(t1)))
            end)
        end)
        nsplit = ceil(Int, (t1 - t0) / timestep(max_gradient))
        append!(control_times, range(t0, t1; length=nsplit+1)[2:end-1])
    end
    sort!(control_times)

    index_ro = ones(length(waveforms))
    instants = map(control_times) do t
        map(1:length(waveforms)) do index_seq
            i_ro = index_ro[index_seq]
            if i_ro <= length(ro) && (t == ro[i_ro].time)
                index_ro[index_seq] += 1
                return ro[i_ro]
            end
            return nothing
        end
    end
    vector_type = static_vector_type(length(waveforms))
    return build_blocks(control_times, vector_type, waveforms, grad_interpolators, instants)
end

function build_blocks(control_times, vector_type, waveforms, grad_interpolators, instants)
    return map(control_times[1:end-1], control_times[2:end], instants[2:end]) do t0, t1, instant
        MultSequencePart{vector_type}(
            t1 - t0,
            map(waveforms, grad_interpolators) do waveform, interpolators
                get_block(waveform, t0, t1, interpolators)
            end,
            instant
        )
    end
end

function get_block(waveform, t0, t1, grad_interpolators)
    @assert t0 <= t1
    if t0 == t1
        return EmptyPart()
    end
    grad0 = SVector{3, Float64}(grad_interpolators[1](t0), grad_interpolators[2](t0), grad_interpolators[3](t0))
    grad1 = SVector{3, Float64}(grad_interpolators[1](t1), grad_interpolators[2](t1), grad_interpolators[3](t1))

    grad_part = if all(iszero.(grad0)) && all(iszero.(grad1))
        EmptyPart()
    elseif all(grad0 .== grad1)
        ConstantPart(grad0) 
    else
        LinearPart(grad0, grad1)
    end

    frac_start = 1.
    frac_final = 1.
    for (t0_rf, t1_rf, rf_pulse) in waveform.rf
        if t0_rf >= t1 || t1_rf <= t0
            continue
        end
        @assert t0_rf <= t0 && t1_rf >= t1 "RF pulse should not start or end within the block. This should have been split in the control times."

        full_pulses = PulsePart(rf_pulse, EmptyPart())
        if t0_rf < t0
            frac_start = (t0 - t0_rf) / (t1_rf - t0_rf)
            _, full_pulses = split_fraction(full_pulses, frac_start)
        end
        if t1_rf > t1
            frac_final = (t1 - t0) / (t1_rf - t0)
            full_pulses, _ = split_fraction(full_pulses, frac_final)
        end
        return PulsePart(full_pulses.pulse, grad_part, full_pulses.first_duration, full_pulses.last_duration)
    end
    return grad_part
end

empty_sequence() = KomaMRIBase.Sequence()


"""
    gradient_waveform(sequence, index::Int)

Extracts the gradient waveform for a sequence.

This needs to be implemented for any sequence type that needs to be processed by the simulator in addition to `get_pulses` and `readout_times`.
"""
function gradient_waveform(sequence::KomaMRIBase.Sequence, index::Int)
    times = [[0.]]
    ampls = [[0.]]
    current_time = 0.
    for block in sequence
        obj = block.GR[index, 1]
        new_times = min.(KomaMRIBase.times(obj), KomaMRIBase.dur(block))
        new_ampls = KomaMRIBase.ampls(obj)
        push!(times, new_times .+ current_time)
        push!(ampls, new_ampls)
        current_time += KomaMRIBase.dur(block)
    end
    push!(times, [KomaMRIBase.dur(sequence)])
    push!(ampls, [0.])
    ftimes = vcat(times...)
    fampls = vcat(ampls...)
    return ftimes, fampls
end


"""
    get_pulses(sequence)

Extracts the pulse events for a sequence.

This needs to be implemented for any sequence type that needs to be processed by the simulator in addition to `gradient_waveform` and `readout_times`.

It returns a tuple with:
- the time of the start of the pulse event in ms since the beginning of the sequence.
- the time of the end of the pulse event in ms since the beginning of the sequence.
- a vector of `ConstantPulse` objects representing the pulse event.
"""
function get_pulses(sequence::KomaMRIBase.Sequence)
    pulses = Tuple{Float64, Float64, Vector{ConstantPulse}}[]
    for (index_block, block) in enumerate(sequence)
        times = KomaMRIBase.times(block.RF[1])
        if length(times) == 0
            continue
        end
        complex_ampls = KomaMRIBase.ampls(block.RF[1])
        ampls = abs.(complex_ampls)
        phases = rad2deg.(angle.(complex_ampls))
        freqs = fill(NaN, length(times))
        for i in 2:length(freqs)-1
            lower_index = i - 1
            upper_index = i + 1
            if iszero(ampls[i])
                continue
            end
            while lower_index > 0 && iszero(ampls[lower_index])
                lower_index -= 1
            end
            if lower_index == 0
                lower_index = i
            end
            while upper_index <= length(ampls) && iszero(ampls[upper_index])
                upper_index += 1
            end
            if upper_index == length(ampls) + 1
                upper_index = i
            end
            if lower_index == upper_index
                error("Frequency can not be determined for isolated non-zero pulse for block $index_block.")
            end
            dphase = (phases[upper_index] - phases[lower_index])
            if dphase > 180
                dphase -= 360
            elseif dphase < -180
                dphase += 360
            end
            freqs[i] = dphase / (times[upper_index] - times[lower_index]) / 360
        end
        first_index = findfirst(!iszero, ampls)
        last_index = findlast(!iszero, ampls)
        all_timestep = diff(times[first_index:last_index])
        timestep = mean(all_timestep)
        @assert all(isapprox.(all_timestep, timestep; rtol=1e-4)) "non-constant timestep in RF pulse for block $index_block."

        pulses_parts = [ConstantPulse(ampls[i], phases[i], freqs[i]) for i in first_index:last_index]
        push!(pulses, (times[first_index] - timestep/2, times[last_index] + timestep/2, pulses_parts))
    end
    return pulses
end


"""
    readout_times(sequence)

Returns the ADC sampling times for the sequence.

This needs to be implemented for any sequence type that needs to be processed by the simulator in addition to `gradient_waveform` and `get_pulses`.
"""
readout_times(sequence::KomaMRIBase.Sequence) = KomaMRIBase.get_adc_sampling_times(sequence)

function compress_timeseries(times::AbstractVector{<:Number}, values::AbstractVector{<:Number}; rtol=1e-4)
    @assert length(times) == length(values)
    ctimes = [times[1]]
    cvalues = [values[1]]

    next_index = 2
    while next_index <= length(times)
        guess_slope = (values[next_index] - cvalues[end]) / (times[next_index] - ctimes[end])
        while (next_index < length(times))  && isapprox(values[next_index + 1] - values[next_index], guess_slope * (times[next_index+1] - times[next_index]); rtol=rtol)
            next_index += 1
        end
        push!(ctimes, times[next_index])
        push!(cvalues, values[next_index])
        next_index += 1
    end
    return ctimes, cvalues
end


"""
    split_fraction(part::SequencePart, fractions::Number(s))

Splits the `SequencePart` into two parts at the given fraction(s) of the duration of the part.

If fraction is a single number a tuple with the two parts is returned.
If fraction is a vector of numbers, a vector of `SequencePart`s is returned of length `length(fractions) + 1`.
The fractions are assumed to be ordered!
"""
split_fraction(::EmptyPart, ::Number) = (EmptyPart(), EmptyPart())
split_fraction(part::ConstantPart, ::Number) = (ConstantPart(part.strength), ConstantPart(part.strength))
split_fraction(part::LinearPart, fraction::Number) = (
    LinearPart(part.start, part.start + fraction * (part.final - part.start)),
    LinearPart(part.start + fraction * (part.final - part.start), part.final)
)




function split_fraction(part::T, fractions::AbstractVector{<:Number}) where {T}
    result = Vector{T}(undef, length(fractions) + 1)
    current_part = part
    current_fraction = 1.
    for (i, fraction) in enumerate(fractions)
        relative_fraction = (fraction - (1 - current_fraction)) / current_fraction
        if relative_fraction < 0 || relative_fraction > 1
            error("Fractions should be between 0 and 1 and ordered.")
        end
        part1, part2 = split_fraction(current_part, relative_fraction)
        result[i] = part1
        current_part = part2
        current_fraction = 1. - fraction
    end
    result[end] = current_part
    return result
end


"""
    split_number(part::SequencePart, n::Int)

Splits the `SequencePart` into `n` equal parts.
"""
function split_number(part, n::Int)
    if n < 1
        error("n should be at least 1.")
    elseif n == 1
        return [part]
    else
        fractions = (1:n-1) ./ n
        return split_fraction(part, fractions)
    end
end


"""
    InstantSequencePart(instants)

A set of `N` instant pulses/gradients/readouts that should be applied to the spins.

Some of the instants might be `nothing`.
"""
struct InstantSequencePart{N, T, ST<:AbstractVector{T}}
    instants :: ST
    function InstantSequencePart{VT}(instants::AbstractVector) where {VT}
        T = Union{typeof.(instants)...}
        N = length(instants)
        return new{N, T, VT{T}}(VT{T}(instants))
    end
end


"""
    MultSequencePart(duration, parts)

A set of N [`SequencePart`](@ref) objects representing overlapping parts of `N` sequences.
"""
struct MultSequencePart{N, T<:SequencePart, ST<:AbstractVector{T}, IT, IST}
    duration :: Float64
    parts :: ST
    instants :: InstantSequencePart{N, IT, IST}
end

function MultSequencePart{VT}(duration::Number, parts::AbstractVector, instants::AbstractVector) where {VT}
    N = length(parts)
    if iszero(N)
        return new{0, EmptyPart, SVector{0, EmptyPart}}(Float64(duration), SVector{0, EmptyPart}())
    end

    base_type = typeof(parts[1])
    if all(p -> p isa base_type, parts)
        T = base_type
    else
        T = SequencePart
    end
    return MultSequencePart{N, T, VT{T}}(
        Float64(duration),
        VT{T}(parts),
        InstantSequencePart{VT}(instants),
    )
end


"""
    get_readouts(sequence[, start_time]; readouts=nothing, nTR=1, skip_TR=0)
    get_readouts(adc_sample_times, TR[, start_time]; readouts=nothing, nTR=1, skip_TR=0)

Returns whether this sequence repeats and when it will be readout.

The return in a tuple with:
1. a boolean indicating whether the sequence is repeating or not, and
2. a vector of `IndexedReadout`s indicating the readouts that will be used in the simulation.

The sequence `adc_sample_times` and `TR` can be given directly, or they can be extracted from the sequence using `sequencee_waveform(sequence, start_time; kwargs...)`.

This can be used to identify which readouts will be (or have been) used in the simulation by running:
`collect(get_readouts(sequence, snapshot.current_time; kwargs...))`
where `snapshot` is the starting snapshot (which has a `current_time` of 0 by default) and `kwargs` are the keyword arguments used in
[`readout`](@ref MCMRSimulator.Evolve.readout) (i.e., `readouts, `nTR`, and `skip_TR`).

By default the readout/ADC objects with the actual sequence definition are used.
These can be overriden by `readouts`, which can be set to a vector of the timings of the readouts within each TR.

# Non-repeating sequences
This is the behaviour if both `nTR` and `skip_TR` are not set by the user.
Any readouts before the `start_time` are ignored.
Any readouts after the `start_time` (whether from the `readouts` keyword or within the sequence definition) are returned.

# Repeating sequences
This is the behaviour if either `nTR` or `skip_TR` or both are set by the user.
If at `start_time` any of the readouts in the current TR have already passed, then the readout will only start in the next TR (unless `skip_TR` is set to -1).
We will skip an additional number of TRs given by `skip_TR` (default: 0).
Then readouts will continue for the number of TRs given by `nTR` (default: 1).
"""
function get_readouts(adc_sample_times::AbstractVector, TR::Number, start_time::Number=0.; readouts=nothing, nTR=nothing, skip_TR=nothing) 
    repeat = !(isnothing(nTR) && isnothing(skip_TR))
    if isnothing(nTR)
        nTR = 1
    end
    if isnothing(skip_TR)
        skip_TR = 0
    end

    if !isnothing(readouts)
        use_readouts = collect(readouts)
        if repeat
            max_ro = maximum(use_readouts)
            if max_ro > TR
                if isapprox(max_ro, TR, rtol=1e-3)
                    @warn "Adjusting maximum readout time ($max_ro) to match the TR of $TR."
                else
                    error("Readouts have been scheduled at $max_ro ms beyond the sequence repetitition time of $TR ms.")
                end
            end
            use_readouts = min.(use_readouts, TR)
        end
    else
        use_readouts = collect(adc_sample_times)
    end

    if !repeat
        return repeat, [IndexedReadout(time, 0, index) for (index, time) in enumerate(use_readouts) if time >= start_time]
    end

    time_in_TR = start_time % TR
    current_TR = floor(Int, start_time / TR) + 1
    first_TR = if time_in_TR > 0 && any(ro -> ro >= time_in_TR, use_readouts)
        current_TR + 1 + skip_TR
    else
        current_TR + skip_TR
    end

    return repeat, [
        IndexedReadout(TR * (tr - 1) + ro, tr, index) 
        for tr in first_TR:first_TR + nTR - 1 
        for (index, ro) in enumerate(use_readouts)
        if (TR * (tr - 1) + ro) >= start_time
    ]
end

function get_readouts(sequence, start_time::Number=0.; kwargs...) 
    obj = sequence_waveform(sequence)
    get_readouts(obj.samples, obj.TR, start_time; kwargs...)
end


end