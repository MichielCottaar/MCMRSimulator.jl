# Introduction
[MCMRSimulator.jl](https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl) allows simulation of MR signal generation using Monte Carlo simulations.
The spin evolution of randomly diffusing particles is tracked under influence of one or more MR sequences.
At present, the simulator allows to model
- Free diffusion and diffusion restricted by [`Walls`](@ref), [`Cylinders`](@ref), [`Spheres`](@ref), and/or a triangular [`Mesh`](@ref).
- R1 and R2 relaxation using global or local R1/R2 parameters
- MR sequences consisting of arbitrary RF pulses, gradients, and readouts
- Off-resonance magnetic field generation by myelinated cylinders or meshes
- Magnetisation transfer between liquid spins and bound spins in membranes.
- Membrane permeability (i.e., exchange)
- Surface relaxation
- Surface tension of membranes causing spins to get temporarily "stuck" when they hit a membrane

Future (potential) features:
- Off-resonance field by iron particles
- Contribution from metabolites (i.e., spectroscopy)

!!! warning "beta"
    This MR simulator is still under very active development and the API might still change substantially at any time!

We use the following units throughout (unless otherwise noted):
- Times are in ms. Equivalently, RF pulse amplitudes and off-resonance magnetic fields are in kHz (i.e., 1/ms).
- Positions are in um. So, gradients are in kHz/um (not mT/m).
- Angles are in degrees (not radians). These are used for phases (of spins and RF pulses) as well as RF pulse flip angles. 
- Susceptibilities are in parts per million (ppm).

## How to get started?
1. If MCMRSimulator is not yet installed, follow the [installation instructions](@ref installation).
2. Look through one of the tutorials. There are two available, depending on which interface you prefer to use:
    - For the command line interface: [CLI tutorial](@ref tutorial_cli)
    - For the julia interface: [Julia tutorial](@ref tutorial_julia)
3. If you want more information on a specific topic, you can check one of the more dedicated sections:
    - [Geometry](@ref geometry)
    - [MRI & collision properties](@ref properties)
    - [Sequences](@ref sequences)
    - Full [API](@ref api)
## Contributors
The original simulator was written by Michiel Cottaar.

Other contributors:
- Zhiyu Zheng
## Movie of spins moving through cylinders
```@raw html
<iframe src="https://ox.cloud.panopto.eu/Panopto/Pages/Embed.aspx?id=b6211751-2743-4bb8-b65a-af5d011a8684&autoplay=true&offerviewer=false&showtitle=false&showbrand=false&captions=false&interactivity=none" style="border: 1px solid #464646;" allowfullscreen allow="autoplay"></iframe>
```
