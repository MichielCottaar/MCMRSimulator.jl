# [Obstructions to free diffusion](@id geometry)
## Defining the geometry
MRSimulator.jl comes with a variety of basic components that can be used to represent various components in the tissue microstructure.

| Component     | Class            | Constructor         |  Dimensionality |
| ------------- | ---------------- | ------------------- |  -------------- |
| infinite walls | [`Wall`](@ref) | [`walls`](@ref)  |  1 |
| hollow infinite cylinder | [`Cylinder`](@ref) |  [`cylinders`](@ref)   |  2 |
| Annulus with inner and outer cylinders | [`Annulus`](@ref) | [`annuli`](@ref)   |  2 |
| Spirals | [`Spiral`](@ref) | [`spirals`](@ref)   |  2 |
| hollow sphere | [`Sphere`](@ref) |  [`spheres`](@ref)   |  3 |
| mesh | [`Mesh`](@ref) | see [Defining a mesh](@ref) |  3 |

The constructors for these components all have a similar interface.
Some expect certain component-specific arguments (e.g., radii for [`spheres`](@ref) and [`cylinders`](@ref).
Some also have component-specific keyword argumetns (e.g., the keywords regarding the off-resonance produces by [Myelinated cylinders](@ref Myelinated_cylinders)).
Finally, they expect a set of keyword arguments that control their location.
These arguments are identicaly across all constructors (although the expected input depends on the dimensionality of the component as listed in the table above):
- `positions`: Set the positions for each generated components
- `repeats`: Set the distance with which all components should be repeated
- `rotation`: Applies a single rotation to the whole system.
Components with a lower dimensionality are defined by default along the x-axis (for dimensionality of 1) or the x-y plane (for dimensionality of 2). 
In other words, the normal of the [`walls`](@ref) point in the x-axis by default, while the [`cylinders`](@ref) point in the z-axis.
Shifts and repeats should only be provided in this lower-dimensional space.
The `rotation` can be used to define these components along other lines/planes.

For example, we can create two base cylinders, which repeeat infinitely by running:
```@example
using MRSimulator
geometry = cylinders(sqrt(0.5), positions=[[0, 0], [1, 1]], repeats=[2, 2])
using CairoMakie # hide
f = plot(PlotPlane(size=4), geometry) # hide
save("regular_cylinders.png", f) # hide
nothing # hide
```  

![Plot showing two cylinders repeating ad infinitum](regular_cylinders.png)

Alternatively, the same configuration could be produced with a single cylinder by providing a `rotation`.
```@example
using MRSimulator
rotation = [
    sqrt(0.5) sqrt(0.5) 0.
    -sqrt(0.5) sqrt(0.5) 0.
    0. 0. 1.
    ]
geometry = cylinders(sqrt(0.5), repeats=[sqrt(2), sqrt(2)], rotation=rotation)
using CairoMakie # hide
f = plot(PlotPlane(size=4), geometry) # hide
save("regular_cylinders2.png", f) # hide
nothing # hide
```  
![Plot showing single cylinders repeating ad infinitum](regular_cylinders2.png)

Myelin can be added to the cylinders, spirals, or annuli as described [here](@ref off_resonance).

A geometry is defined by either the [`TransformObstruction`](@ref) returned by a single call to these constructors
or by an array of [`TransformObstruction`](@ref) objects.
### Randomly distributed cylinders/annuli/spirals
A random set of positions and radii can be created using [`random_positions_radii`](@ref).
The user in this case sets a target density (70% in the example below) and over which length scale the configuration should repeat itself (20x20 micrometer in the example below).
```@example random_distribution
using MRSimulator # hide
(positions, outer_radii) = random_positions_radii((20, 20), 0.7, 2)
nothing # hide
```

These can be used to produce randomly distributed cylinders:
```@example random_distribution
geometry = cylinders(outer_radii; positions=positions, repeats=(20, 20))
using CairoMakie # hide
f = plot(PlotPlane(size=20), geometry) # hide
save("random_cylinders.png", f) # hide
nothing # hide
```
![Illustrating configuration of random cylinders](random_cylinders.png)

When used as initialisation for annuli or spirals, an inner radius will also need to be computed:
```@example random_distribution
inner_radii = 0.8 .* outer_radii
geometry = annuli(inner_radii, outer_radii; positions=positions, repeats=(20, 20))
using CairoMakie # hide
f = plot(PlotPlane(size=20), geometry) # hide
save("random_annuli.png", f) # hide
nothing # hide
```
![Illustrating configuration of random annuli](random_annuli.png)


## Defining a mesh