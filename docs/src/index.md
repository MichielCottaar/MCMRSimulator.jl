# Introduction
[MRSimulator.jl](https://git.fmrib.ox.ac.uk/ndcn0236/MRSimulator.jl) allows simulation of MR signal generation using Monte Carlo simulations.
The spin evolution of randomly diffusing particles is tracked under influence of one or more MR sequences.
At present, the simulator allows to model
- Free diffusion and diffusion restricted by walls, cylinders, spheres, or meshes
- T1 and T2 relaxation
- MR sequences consisting of RF pulses, gradients, and readouts
- [Off-resonance field](@ref) generation by myelinated cylinders

```@index
```
