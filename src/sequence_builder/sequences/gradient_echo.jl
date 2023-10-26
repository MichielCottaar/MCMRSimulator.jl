module GradientEcho
import ...DefineSequence: define_sequence
import ....Scanners: Scanner
import ....Sequences: InstantRFPulse, Readout, RFPulse
import ...BuildingBlocks: duration

"""
    gradient_echo(TE; TR=<TE>, scanner=<3T scanner>, excitation_pulse=Instant, excitation_time=<half pulse duration>, readout_time=0, crusher=0.)

Creates a gradient echo sequence consisting of:
- `excitation_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `excitation_time` should be set as well.
- a delay
- a readout `TE` ms after the excitation.
- a delay of readout_time/2
- crusher gradient (default: instant with q-value of 1 rad/um). Set `crusher` to a positive number to change the duration or to any MR gradient (see [`gen_crusher`](@ref)).
"""
function gradient_echo(
    TE;
    TR=nothing,
    scanner=Scanner(B0=3.),
    excitation_pulse=InstantRFPulse(flip_angle=90),
    excitation_time=nothing,
    readout_time=0.,
    crusher=0.,
)
    excitation_time = isnothing(excitation_time) ? duration(excitation_pulse) / 2 : excitation_time
    if crusher isa Number
        crusher = gen_crusher(duration=crusher, scanner=scanner)
    end
    define_sequence(scanner, TR) do 
        [
            excitation_pulse,
            TE - duration(excitation_pulse) + excitation_time,
            Readout(),
            max(readout_time/2, 1e-6),
            crusher,
            1e-6
        ]
    end
end
end