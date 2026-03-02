"""
Module defining the building blocks of a sequence, including gradient waveforms, RF pulses and ADC events.
"""
module Base

import StaticArrays: SVector
import ...Scanners: Scanner, B0


"""
Gradient waveform, containing the time points and the corresponding gradient amplitudes in 3D.
"""
struct GradientWaveform
    time::Vector{Float64}
    amplitude::Vector{SVector{3, Float64}}
end

duration(gradient::GradientWaveform) = gradient.time[end]


"""
RF pulse, containing the time points, the corresponding RF amplitudes and phases.
"""
struct RFPulse
    time::Vector{Float64}
    amplitude::Vector{Float64}
    phase::Vector{Float64}
end

duration(rf_pulse::RFPulse) = rf_pulse.time[end]


"""
ADC event, containing the time points of the ADC samples.
"""
struct ADC
    samples::Vector{Float64}
end

duration(adc::ADC) = adc.samples[end]

duration(::Nothing) = 0.


"""
Basic building block of a sequence, containing the duration of the block, the gradient waveform, the RF pulse and the ADC events.
"""
struct BuildingBlock
    duration::Float64
    gradient::Union{GradientWaveform, Nothing}
    rf_pulse::Union{RFPulse, Nothing}
    adc::Union{ADC, Nothing}
end

duration(block::BuildingBlock) = block.duration


"""
Sequence, containing the building blocks of a sequence, including gradient waveforms, RF pulses and ADC events.
"""
struct Sequence
    blocks::Vector{BuildingBlock}
    scanner::Scanner
end

B0(sequence::Sequence) = B0(sequence.scanner)


"""
    duration(sequence/block)

Returns the duration of a sequence or a building block in milliseconds.
"""
duration(sequence::Sequence) = sum(duration.(sequence.blocks))

end