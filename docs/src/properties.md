# [Simulation properties](@id properties)
How the spins behave is determined by the tissue [geometry](@ref geometry), the applied [MRI sequences](@ref sequence), and user-provided flags determining how the spin magnetisation evolves. Here we discuss how the spin magnetisation evolution can be affected by these user-provided flags.

For example, one such flag is the `diffusivity`, which can be set as a keyword argument while generating the [`Simulation`](@ref).
## MRI properties
[`MRIProperties`](@ref) determine the spin evolution for free and stuck particles. They include:
- the longitudinal relaxation time `T1` or rate `R1`
- the transverse relaxation time `T2` or rate `R2`
- the global `off_resonance` field (i.e., any off-resonance not caused by the sequence or the geometry)
At the [`Simulation`](@ref) level these parameters can be set by supplying the `T1`/`R1`, `T2`/`R2`, or `off_resonance` flags (see [`GlobalProperties`](@ref)), such as:
```julia
simulation = Simulation(sequences, T1=2000, R2=1/80)
```
These MRI properties can be locally overwritten when defining the [geometry](@ref geometry) (see [`ObstructionProperties`](@ref)). The relaxation properties within a compartment can be set by setting `T1_inside`/`R1_inside`, `T2_inside`/`R2_inside`, and/or `off_resonance_inside` flags. Equivalently, the relaxation properties for particles stuck on the surface can be set using the `T1_surface`/`R1_surface`, `T2_surface`/`R2_surface`, and/or `off_resonance_surface` flags.

If not set at the global or local level, there will be no longitudinal or transverse relaxation and there will be no off-resonance field.

## Collision properties
[`CollisionProperties`](@ref) determine the behaviour of spins at the time of a collision. Like MRI properties they can be set at the global level ([`GlobalProperties`](@ref)) or overwritten at the local level ([`ObstructionProperties`](@ref)). There are four such properties:
- [`MT_fraction`](@ref): the fraction of transverse signal lost (and longitudinal signal regained) at every collision. This fraction is adjusted to take into account the timestep (see [`correct_for_timestep`](@ref)). Note that this is not the recommended way to model magnetisation transfer. Instead, we recommend using the `surface_density` as discussed below.
- [`permeability`](@ref): the probability of the spin passing through the surface. If the spins do not pass through, they will undergo regular reflection (or get stuck, see below). Like `MT_fraction` it will be adjusted to take into account the timestep (see [`correct_for_timestep`](@ref)).
- [`surface_density`](@ref) and [`dwell_time`](@ref): These control the density and dwell time of spins on the surface. Depending on the MRI properties assigned to these stuck particles (see above), these stuck particles can be used to represent water stuck at the membranes due to surface tension or spins in the membrane itself (which is in exchange with the free water through magnetisation transfer).

