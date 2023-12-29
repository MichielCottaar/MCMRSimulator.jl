module StimulatedEcho
import ....Scanners: Scanner, max_gradient
import ....Sequences: InstantRFPulse, Readout, RFPulse, InstantGradient, MRGradients, rotate_bvec
import ...DefineSequence: define_sequence
import ...Diffusion: add_linear_diffusion_weighting
import ...BuildingBlocks: duration
import StaticArrays: SVector

"""
stimulated_echo(TE, 
    stimulate_interval;
    TR=nothing,
    scanner=Scanner(B0=3.),
    excitation_pulse=InstantRFPulse(flip_angle=90),
    excitation_time=nothing,
    refocus_pulse=InstantRFPulse(flip_angle=180, phase=90),
    refocus_time=nothing,
)

Creates a stimulated echo sequence consisting of:
- `excitation_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `excitation_time` should be set as well.
- a delay
- `stimulate_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `stimulate_time` should be set as well.
- a readout `TE` ms after the excitation.

"""

function stimulated_echo(TE, 
    stimulate_interval;
    TR=nothing,
    scanner=Scanner(B0=3.),
    excitation_pulse=InstantRFPulse(flip_angle=90),
    excitation_time=nothing,
    stimulate_pulse=InstantRFPulse(flip_angle=90, phase=90),
    stimulate_time=nothing,
    spoiler_orientation=SVector{3, Float64}([1., 0., 0.])
)
    if isinf(max_gradient(scanner)) || isinf(max_slew_rate(scanner)) # Maybe implement a check for only inf gradient max in scanner()?
        ramp_time = 0.0000001 # Non-zero to avoid timing issues
    else
        ramp_time = max_gradient(scanner) / max_slew_rate(scanner)
    end

    interval_spoiler = rotate_bvec(MRGradients([
        (0, 0.), 
        (ramp_time, max_gradient(scanner)),
        (stimulate_interval - ramp_time - duration(stimulate_pulse), max_gradient(scanner)),
        (stimulate_interval - duration(stimulate_pulse), 0.), 
    ], apply_bvec=true), spoiler_orientation)

    excitation_time = isnothing(excitation_time) ? duration(excitation_pulse) / 2 : excitation_time
    stimulate_time = isnothing(stimulate_time) ? duration(stimulate_pulse) / 2 : stimulate_time
    define_sequence(scanner, TR) do 
        [
            excitation_pulse,
            TE/2 - duration(excitation_pulse) + excitation_time - stimulate_time,
            stimulate_pulse,
            interval_spoiler,
            stimulate_pulse,
            TE/2 - duration(stimulate_pulse) + stimulate_time,
            Readout()
        ]
    end
end

"""
    dwste(stimulate_interval; 
        TE=20., TR=<TE>+<stimulate_interval>, scanner=<3T scanner>, 
        excitation_pulse=Instant, excitation_time=<half pulse duration>, 
        stimulate_pulse=Instant, stimulate_time=<half pulse duration>,
        readout_time=0, diffusion_time=<maximum>, gradient_duration=<maximum>,
        bval/qval/gradient_strength=<one is required>, orientation=:x,
        )

Creates a stimulated echo sequence consisting of:
- `excitation_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `excitation_time` should be set as well.
- diffusion weighting
- 2x `stimulate_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `refocus_time` should be set as well.
- 'stimulate_interval' between 2 stimulate_pulses
- identical diffusion weighting cancelling out the first one
- a readout `TE` ms after the excitation.
The refocus time is always halfway the excitation time and the readout.

The timings of the RF pulses is set by `TE` and `TR`. The gradient timings will also be affected by `gradient_duration`, `diffusion_time`, `scanner` (which sets the ramp time) and readout_time:
- By default the gradient durations are set to the maximum value possible within the echo time (`TE`) keeping in mind the time needed for the MR readout (`readout_time`) and the time needed to ramp to the maximum gradient strength (set by the `scanner`).
- When gradient_duration is set to 0, the gradient pulses are assumed to be instanteneous (i.e., using [`InstantGradient`](@ref)). The time between these instant gradients can be set using diffusion_time (defaults to TE/2).

The strength of the diffusion gradients is set by one of `bval` (units: ms/um^2), `qval` (units: 1/um), or `gradient_strength` (units: kHz/um). If the resulting gradient strength exceeds the maximum allowed for the `scanner` an AssertionError is raised. The gradient orientation is set by `orientation`.

For more details of how the diffusion weighting is inserted in the [`spin_echo`](@ref) sequence, see [`add_linear_diffusion_weighting`](@ref).
"""
function dwste(stimulate_interval;
    TE=20.,
    TR=nothing,
    scanner=Scanner(B0=3.),
    excitation_pulse=InstantRFPulse(flip_angle=90, phase=90),
    excitation_time=nothing,
    stimulate_pulse=InstantRFPulse(flip_angle=90, phase=90),
    stimulate_time=nothing,
    bval=nothing,
    diffusion_time=nothing,
    qval=nothing,
    gradient_strength=nothing, # in mT/m
    gradient_duration=nothing,
    readout_time=0.,
    orientation=SVector{3, Float64}([1., 0., 0.]),
)

    define_sequence(scanner, TR) do
        sequence = stimulated_echo(
            TE,
            stimulate_interval;
            excitation_pulse=excitation_pulse, excitation_time=excitation_time,
            stimulate_pulse=stimulate_pulse, stimulate_time=stimulate_time, scanner=scanner, spoiler_orientation=orientation
        )

        if !iszero(readout_time)
            sequence[6] -= readout_time/2
            sequence[7] = [
                readout_time/2,
                Readout(),
                readout_time/2
            ]
            sequence[8] -= readout_time/2
        end

        return add_linear_diffusion_weighting(
            sequence, 2, 6,
            bval=bval, diffusion_time=diffusion_time, qval=qval,
            gradient_strength=gradient_strength, gradient_duration=gradient_duration,
            orientation=orientation, scanner=scanner
        )
    end
end

end
