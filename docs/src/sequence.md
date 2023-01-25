# [MR sequences](@id sequence)
## Built-in MR sequences
### Diffusion-weighted MRI
A pulsed-gradient spin-echo can be created using [`dwi`](@ref)
```@example
using MCMRSimulator # hide
sequence = dwi(TR=100., bval=3.)
using CairoMakie # hide
f = plot(sequence) # hide
save("dwi_sequence.png", f) # hide
nothing # hide
```  
![DWI sequence diagram](dwi_sequence.png)

```@example
using MCMRSimulator # hide
sequence = dwi(TR=100., bval=3., gradient_duration=0.)
using CairoMakie # hide
f = plot(sequence) # hide
save("dwi_sequence_delta_0.png", f) # hide
nothing # hide
```  
![DWI sequence diagram with instanteneous gradients](dwi_sequence_delta_0.png)
## Custom MR sequences
In MCMRSimulator.jl an MR [`Sequence`](@ref) describes the RF pulses and gradients applied by the MRI scanner.

This sequence contains:
- A set of instantaneous RF pulses ([`InstantRFPulse`](@ref)), gradients ([`InstantGradient`](@ref)), or readouts ([`Readout`](@ref)). Each of these changes or reads the spin orientations at a specific timepoint.
- A gradient profile from [`create_gradients`](@ref). They can be rotated using [`rotate_bvec`](@ref).
- TODO: how the RF pulse changes over time
Each of these sequence components will play identically every repetition time (TR) of the sequence.

## Defining the MR gradients
The MRI scanner gradients cause the spins to precess at different rates in different part of the tissue.
This encodes the spin location into the spin orientation and can hence be used to measure the movement of spins (e.g., in diffusion-weighted MRI).
Crusher gradients are also commonly used to get rid of unwanted signal contributions due to imperfect RF pulses.

The gradient profile over time can be modeled as:
- Stepwise function: [`StepWiseGradients`](@ref).
- Gradient profile consisting only of straight lines: [`LinearGradients`](@ref).
- TODO: Generic gradient profiles
All of these can be created using [`create_gradients`](@ref).
The off-resonance field at a specific time (or averaged over a timespan) and position can be computed using [`get_gradient`](@ref).
