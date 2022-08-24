# Obstructions to free diffusion
## Defining the geometry
MRSimulator.jl comes with a variety of basic components that can be used to represent various components in the tissue microstructure.

| Component     | Class            | Constructor         | dimensionality |
| ------------- | ---------------- | ------------------- | -------------- |
| infinite walls | [`Wall`](@ref) | [`walls`](@ref)   | 1 |
| hollow infinite cylinder | [`Cylinder`](@ref) | [`cylinders`](@ref)   | 2 |
| Annulus with inner and outer radius | [`Annulus`](@ref) | [`annuli`](@ref)   | 2 |
| hollow sphere | [`Sphere`](@ref) | [`spheres`](@ref)   | 3 |
| mesh | [`Mesh`](@ref) | see [Defining a mesh](@ref)   | 3 |

The constructors for these components all have a similar interface.
Some expect certain component-specific arguments (e.g., radii for [`spheres`](@ref) and [`cylinders`](@ref).
Some also have component-specific keyword argumetns (e.g., the keywords regarding the off-resonance produces by [Myelinated cylinders](@ref)).
Finally, they expect a set of keyword arguments that control their location.
These arguments are identicaly across all constructors (although the expected input depends on the dimensionality of the component as listed in the table above):
- `positions`: Set the positions for each generated components
- `repeats`: Set the distance with which all components should be repeated
- `rotation`: Applies a single rotation to the whole system. This can be an object from 
Components with a lower dimensionality are defined by default along the x-axis (for dimensionality of 1) or the x-y plane (for dimensionality of 2). 
In other words, the normal of the [`walls`](@ref) point in the x-axis by default, while the [`cylinders`](@ref) point in the z-axis.
Shifts and repeats should only be provided in this lower-dimensional space.
The `rotation` can be used to define these components along other lines/planes.

For example, we can create two base spheres, which repeeat infinitely by running:
```julia
import MRSimulator as mr
geometry = mr.spheres(sqrt(3), positions=[[0, 0, 0], [1, 1, 1]], repeats=[2, 2, 2])
```  
This in effect creates a tightly packed, infinitely repeating crystal of hollow spheres in a body-centered cubic arrangement.

A geometry can be formed by one or more calls to these constructors.
This can be passed on to [`Simulation`](@ref).


## Defining a mesh