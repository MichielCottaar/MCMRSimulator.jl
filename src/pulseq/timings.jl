"""
Extract timings, amplitudes, and other information from a Pulseq sequence.
"""
module Timings
import ..Types: PulseqSequence, PulseqBlock, PulseqADC, PulseqGradient, PulseqTrapezoid, PulseqRFPulse
import Interpolations: linear_interpolation
import StaticArrays: SVector


"""
    duration(sequence, unit=:second(default)/:ms/:raster)
    duration(block)

Returns the duration of a sequence in seconds (default), milliseconds, or in units of the block duration raster.

The duration of a building block is always returned in units of the block duration raster.
"""
function duration(sequence::PulseqSequence, unit=:second)
    mult_unit = if unit == :raster
        1
    elseif unit in (:s, :second)
        sequence.definitions.BlockDurationRaster
    elseif unit in (:ms, :millisecond)
        sequence.definitions.BlockDurationRaster * 1e3
    else
        error("Unknown unit: $unit. Use :second, :ms, or :raster.")
    end
    return sum(duration, sequence.blocks) * mult_unit
end

duration(block::PulseqBlock) = block.duration


"""
    adc_sample_times(sequence, unit=:second)
    adc_sample_times(block, adc_raster_time)
    adc_sample_times(adc_object, adc_raster_time)

Returns the times at which ADC samples are taken in a sequence, building block, or ADC object.

If `unit` is set to :raster, the times are returned as tuples of (block_time, adc_time) in units of the block duration raster and seconds, respectively. This is the only supported bevahiour for blocks and ADC objects.
If `unit` is set to :second (default for sequences) or :ms, the times are returned in seconds or milliseconds, respectively.
"""
function adc_sample_times(sequence::PulseqSequence, unit=:second)
    if unit == :raster
        times = Tuple{Int, Float64}[]
        current_time = 0
        for block in sequence.blocks
            for adc_time in adc_sample_times(block, sequence.definitions.AdcRasterTime)
                push!(times, (current_time, adc_time))
            end
            current_time += duration(block)
        end
        return times
    elseif unit in (:s, :second)
        base = adc_sample_times(sequence, :raster)
        return map(base) do (block_time, adc_time)
            block_time * sequence.definitions.BlockDurationRaster + adc_time
        end
    elseif unit in (:ms, :millisecond)
        adc_sample_times(sequence, :second) .* 1e3
    else
        error("Unknown unit: $unit. Use :second, :ms, or :raster.")
    end 
end

adc_sample_times(block::PulseqBlock, adc_raster::Number) = adc_sample_times(block.adc, adc_raster)

function adc_sample_times(adc::PulseqADC, adc_raster::Number)
    edge_times = @. 1e-6 * adc.delay + 1e-9 .* (adc.dwell * (0:adc.num))
    edge_on_raster = round.(edge_times / adc_raster) * adc_raster
    if !all(isapprox.(edge_times, edge_on_raster; atol=adc_raster/100))
        @warn "ADC sample times are not exactly on the ADC raster. Consider adjusting the ADC delay or dwell time to align with the raster."
    end
    return (edge_on_raster[1:end-1] + edge_on_raster[2:end]) / 2
end
adc_sample_times(::Nothing, ::Number) = Float64[]


"""
    gradient_waveform(sequence, dimension=:all, time_unit=:second, amplitude_unit=:Hz_per_m)
    gradient_waveform(block, dimension=:all, gradient_raster_time)
    gradient_waveform(pulseq_gradient, gradient_raster_time)

Returns the gradient waveform of a sequence, building block, or gradient object.

The `dimension` can be set to :x, :y, :z, or :all (default) to return the gradient waveform for the corresponding dimension(s).
The `time_unit` can be set to :second (default) or :ms to return the times in seconds or milliseconds, respectively. This is only supported for sequences.
For blocks and gradient objects, the times are always returned in units of the block duration raster/2.

Gradient amplitudes are returned in units of Hz/m. They can be converted to `kHz/um` (as used in the simulator) by multiplying by 1e3.
"""
gradient_waveform(object, dimension::Symbol, args...) = gradient_waveform(object, Val(dimension), args...)
gradient_waveform(object, dimension::Int, args...) = gradient_waveform(object, Val((:x, :y, :z)[dimension]), args...)

function gradient_waveform(grad::PulseqBlock, ::Val{:all}, gradient_raster_time::Number)
    t1, g1 = gradient_waveform(grad, Val(:x), gradient_raster_time)
    t2, g2 = gradient_waveform(grad, Val(:y), gradient_raster_time)
    t3, g3 = gradient_waveform(grad, Val(:z), gradient_raster_time)
    all_times = sort!(unique(vcat(t1, t2, t3)))
    if length(all_times) == 0
        return all_times, SVector{3, Float64}[]
    end
    interpolators = (
        length(t1) == 0 ? t -> 0. : linear_interpolation(t1, g1, extrapolation_bc=0),
        length(t2) == 0 ? t -> 0. : linear_interpolation(t2, g2, extrapolation_bc=0),
        length(t3) == 0 ? t -> 0. : linear_interpolation(t3, g3, extrapolation_bc=0),
    )
    return all_times, [SVector{3, Float64}(interp(time) for interp in interpolators) for time in all_times]
end

for symb in (:x, :y, :z)
    sgx = Symbol("g", symb)
    @eval function gradient_waveform(grad::PulseqBlock, ::Val{$(QuoteNode(symb))}, gradient_raster_time::Number)
        gradient_waveform(grad.$sgx, gradient_raster_time)
    end
end

function gradient_waveform(grad::PulseqTrapezoid, gradient_raster_time::Number)
    ndelay = round(Int, grad.delay * 1e-6 / gradient_raster_time)
    @assert ndelay ≈ grad.delay * 1e-6 / gradient_raster_time "Trapezoid delay is not an integer multiple of the gradient raster time."
    durations = [ndelay, grad.rise, grad.flat, grad.fall]
    times = cumsum(durations)
    amps = [0., grad.amplitude, grad.amplitude, 0.]
    return times, amps
end

function gradient_waveform(grad::PulseqGradient, gradient_raster_time::Number)
    ndelay = round(Int, grad.delay * 1e-6 / gradient_raster_time)
    @assert ndelay ≈ grad.delay * 1e-6 / gradient_raster_time "Gradient delay is not an integer multiple of the gradient raster time."
    if isnothing(grad.time)
        times = collect(1:2:(2 * length(grad.shape)))
        ampls = grad.shape.samples .* grad.amplitude
        prepend!(times, 0)
        prepend!(ampls, 0.)
        append!(times, 2 * length(grad.shape))
        append!(ampls, 0.)
        return ndelay .+ times, ampls
    else
        times = grad.time * 2
        ampls = grad.shape.samples .* grad.amplitude
        return ndelay .+ times, ampls
    end
end

gradient_waveform(::Nothing) = (Int[], Float64[])

gradient_waveform(seq::PulseqSequence) = gradient_waveform(seq, Val(:all), :second)

function gradient_waveform(seq::PulseqSequence, dimension::Val{D}, time_unit::Symbol=:second) where {D}
    if time_unit == :raster
        times = Tuple{Int, Int}[]
        amplitude = (D == :all ? SVector{3, Float64} : Float64)[]
        current_time = 0
        for block in seq.blocks
            (block_time, grad) = gradient_waveform(block, dimension)
            for sub_time in block_time
                push!(times, (current_time, sub_time))
            end
            append!(amplitude, grad)
            current_time += duration(block)
        end
        return times, amplitude
    elseif time_unit in (:s, :second)
        base, grads = gradient_waveform(seq, dimension, :raster)
        times = map(base) do (block_time, sub_time)
            block_time * seq.definitions.BlockDurationRaster + sub_time * seq.definitions.GradientRasterTime / 2
        end
        return times, grads
    elseif time_unit in (:ms, :millisecond)
        times, grads = gradient_waveform(seq, dimension, :second)
        return times .* 1e3, grads
    else
        error("Unknown time unit: $time_unit. Use :second, :ms, or :raster.")
    end
end


"""
    rf_pulses(sequence, unit=:second)

Returns a list of all RF pulses in a sequence.

Each pulse is returned as a tuple of the:
- The pulse start time (including the pulse delay). If `unit` is set to :second (default) or :ms, the time is returned in seconds or milliseconds, respectively.
- `PulseqRFPulse` object.
"""
function rf_pulses(seq::PulseqSequence, unit=:second)
    pulses = Tuple{unit == :raster ? Int : Float64, PulseqRFPulse}[]
    current_time = 0
    for block in seq.blocks
        if !isnothing(block.rf)
            pulse = block.rf
            if unit in (:s, :second)
                push!(pulses, (current_time * seq.definitions.BlockDurationRaster + pulse.delay * 1e-6, pulse))
            elseif unit in (:ms, :millisecond)
                push!(pulses, ((current_time * seq.definitions.BlockDurationRaster + pulse.delay * 1e-6) * 1e3, pulse))
            else
                error("Unknown unit: $unit. Use :second, :ms, or :raster.")
            end
        end
        current_time += duration(block)
    end
    return sort(pulses; by=first)
end



end