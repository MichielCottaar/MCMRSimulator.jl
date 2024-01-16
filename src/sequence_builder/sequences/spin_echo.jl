module SpinEcho
import ....Scanners: Scanner, max_gradient
import ....Sequences: InstantRFPulse, Readout, RFPulse, InstantGradient, rotate_bvec
import ...DefineSequence: define_sequence
import ...Diffusion: add_linear_diffusion_weighting, gen_crusher
import ...BuildingBlocks: duration
import StaticArrays: SVector
"""
    spin_echo(TE; TR=<TE>, scanner=<3T scanner>, excitation_pulse=Instant, excitation_time=<half pulse duration>, refocus_pulse=Instant, refocus_time=<half pulse duration>)

Creates a gradient echo sequence consisting of:
- `excitation_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `excitation_time` should be set as well.
- a delay
- `refocus_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the refocus does not take place halfway the RF pulse, `refocus_time` should be set as well.
- a readout `TE` ms after the excitation.
- a delay of readout_time/2
- crusher gradient (default: instant with q-value of 1 rad/um). Set `crusher` to a positive number to change the duration or to any MR gradient (see [`gen_crusher`](@ref)).
The refocus time is always halfway the excitation time and the readout.
"""
function spin_echo(TE;
    TR=nothing,
    scanner=Scanner(B0=3.),
    excitation_pulse=InstantRFPulse(flip_angle=90),
    excitation_time=nothing,
    refocus_pulse=InstantRFPulse(flip_angle=180, phase=90),
    refocus_time=nothing,
    readout_time=0.,
    crusher=0.,
)
    excitation_time = isnothing(excitation_time) ? duration(excitation_pulse) / 2 : excitation_time
    refocus_time = isnothing(refocus_time) ? duration(refocus_pulse) / 2 : refocus_time
    if crusher isa Number
        crusher = gen_crusher(duration=crusher, scanner=scanner)
    end
    define_sequence(scanner, TR) do 
        [
            excitation_pulse,
            TE/2 - duration(excitation_pulse) + excitation_time - refocus_time,
            refocus_pulse,
            TE/2 - duration(refocus_pulse) + refocus_time - readout_time/2,
            readout_time/2,
            Readout(),
            max(readout_time/2, 1e-6),
            crusher,
            1e-6
        ]
    end
end

"""
    dwi(; 
        TE=80., TR=<TE>, scanner=<3T scanner>, 
        excitation_pulse=Instant, excitation_time=<half pulse duration>, 
        refocus_pulse=Instant, refocus_time=<half pulse duration>,
        readout_time=0, diffusion_time=<maximum>, gradient_duration=<maximum>,
        bval/qval/gradient_strength=<one is required>, orientation=:x,
        )

Creates a gradient echo sequence consisting of:
- `excitation_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the excitation does not take place halfway the RF pulse, `excitation_time` should be set as well.
- diffusion weighting
- `refocus_pulse`: by default this is an [`InstantRFPulse`](@ref), but can be replaced with an [`RFPulse`](@ref). If the refocus does not take place halfway the RF pulse, `refocus_time` should be set as well.
- identical diffusion weighting cancelling out the first one
- a readout `TE` ms after the excitation.
The refocus time is always halfway the excitation time and the readout.

The timings of the RF pulses is set by `TE` and `TR`. The gradient timings will also be affected by `gradient_duration`, `diffusion_time`, `scanner` (which sets the ramp time) and readout_time:
- By default the gradient durations are set to the maximum value possible within the echo time (`TE`) keeping in mind the time needed for the MR readout (`readout_time`) and the time needed to ramp to the maximum gradient strength (set by the `scanner`).
- When gradient_duration is set to 0, the gradient pulses are assumed to be instanteneous (i.e., using [`InstantGradient`](@ref)). The time between these instant gradients can be set using diffusion_time (defaults to TE/2).

The strength of the diffusion gradients is set by one of `bval` (units: ms/um^2), `qval` (units: 1/um), or `gradient_strength` (units: kHz/um). If the resulting gradient strength exceeds the maximum allowed for the `scanner` an AssertionError is raised. 

The gradient orientation is set by `orientation`.  If the gradient orientation is not set during construction, it can be later applied using [`rotate_bvec`](@ref).

For more details of how the diffusion weighting is inserted in the [`spin_echo`](@ref) sequence, see [`add_linear_diffusion_weighting`](@ref).
"""
function dwi(;
    TE=80.,
    TR=nothing,
    scanner=Scanner(B0=3.),
    excitation_pulse=InstantRFPulse(flip_angle=90, phase=90),
    excitation_time=nothing,
    refocus_pulse=InstantRFPulse(flip_angle=180, phase=0),
    refocus_time=nothing,
    bval=nothing,
    diffusion_time=nothing,
    qval=nothing,
    gradient_strength=nothing, # in mT/m
    gradient_duration=nothing,
    readout_time=0.,
    orientation=:x,
    crusher=0.
)
    define_sequence(scanner, TR) do
        sequence = spin_echo(
            TE;
            excitation_pulse=excitation_pulse, excitation_time=excitation_time,
            refocus_pulse=refocus_pulse, refocus_time=refocus_time, scanner=scanner,
            readout_time=readout_time, crusher=crusher
        )

        return add_linear_diffusion_weighting(
            sequence, 2, 4,
            bval=bval, diffusion_time=diffusion_time, qval=qval,
            gradient_strength=gradient_strength, gradient_duration=gradient_duration,
            orientation=orientation, scanner=scanner
        )
    end
end

end