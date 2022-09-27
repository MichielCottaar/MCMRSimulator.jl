# MR sequences
## Built-in MR sequences
- Diffusion-weighted MRI: [dwi](@ref)
## Custom MR sequences
In MRSimulator.jl an MR [`Sequence`](@ref) describes the RF pulses and gradients applied by the MRI scanner.

This sequence contains:
- A set of instantaneous RF pulses ([`RFPulse`](@ref)), gradients ([`InstantGradient`](@ref)), or readouts ([`Readout`](@ref)). Each of these changes or reads the spin orientations at a specific timepoint.
- A gradient profile from [`create_gradients`](@ref). They can be rotated using [`rotate_bvec`](@ref).
- TODO: how the RF pulse changes over time
Each of these sequence components will play identically every repetition time (TR) of the sequence.

## Defining the MR gradients
The MRI scanner gradients cause the spins to precess at different rates in different part of the tissue.
This encodes the spin location into the spin orientation and can hence be used to measure the movement of spins (e.g., in diffusion-weighted MRI).
Crusher gradients are also commonly used to get rid of unwanted signal contributions due to imperfect RF pulses.

The gradient profile over time can be modeled as:
- Stepwise function: [`StepGradients`](@ref).
- Gradient profile consisting only of straight lines: [`LinearGradients`](@ref).
- Generic gradient profiles: [`GenericGradients`](@ref).
All of these can be created using [`create_gradients`](@ref).
The off-resonance field at a specific time (or averaged over a timespan) and position can be computed using [`get_gradient`](@ref).
