# [Liminal geometry](@id liminal_geometry)

Liminal geometry models a population of cells without explicitly constructing a globally packed geometry. It is useful when the MRI signal depends on the morphology and orientation of individual cells, but a physically packed arrangement of all cells would be expensive or impossible to generate.

## Why liminal geometry?

Explicitly packing cells requires solving a global geometry problem: every cell must be positioned while avoiding unwanted overlaps and preserving the desired volume fractions. This problem becomes particularly difficult for complex morphologies, such as branching cells loaded from SWC files. It can also become the dominant computational cost before any MRI simulation has started.

Liminal geometry avoids this global packing problem. It describes the tissue as a statistical population of independently sampled cell instances. Each instance uses one of the supplied cell geometries, and the simulator applies random translations when it needs to represent an encounter with a cell. The MRI signal can therefore reflect arbitrary cell morphologies and compartment populations without first finding one globally consistent packing of those cells.

This is the central benefit of the model: it separates the local geometry needed to calculate encounters from the global packing problem that is often difficult to solve. The approach is especially useful for heterogeneous populations, mixtures of cell morphologies, and morphologies obtained directly from experimental data.

## Implementation

A liminal geometry contains two types of space:

- an extracellular liminal space, represented statistically rather than by an explicit boundary;
- intracellular regions generated from the child geometries supplied by the user.

The `extracellular_fraction` specifies the volume fraction assigned to the extracellular space. The remaining volume is assigned to the intracellular cell population. When a spin is in the extracellular space, encounters with cells are generated from the population's surface statistics. Once a spin enters a cell, the child geometry controls its local collisions and relaxation properties.

The mean free path is generally direction-dependent. An elongated cell, for example, has a different projected surface area for motion parallel to its axis than for motion perpendicular to it. Consequently, the encounter rate and mean free path can vary with direction. This gives liminal geometry an orientation-dependent effective tortuosity, even though no globally packed arrangement is explicitly stored.

## Julia usage

Create a liminal geometry by supplying a list of `(number_fraction, geometry)` pairs:

```julia
using MCMRSimulator

geometry = LiminalGeometry(
    geometries=[
        (0.7, Spheres(radius=1.0)),
        (0.3, Cylinders(radius=0.8, rotation=:x)),
    ],
    extracellular_fraction=0.2,
)

fixed_geometry = fix(geometry)
```

The child geometries can be any supported geometry that can be used as a three-dimensional cell template, including morphologies loaded from files:

```julia
cell_morphology = read_geometry("cell.swc")
geometry = LiminalGeometry(
    geometries=[(1.0, cell_morphology)],
    extracellular_fraction=0.2,
)
```

The number fractions describe the relative abundance of the cell types. They are not volume fractions. For example, two cell types with number fractions `0.5` and `0.5` are equally common, even if one cell type has twice the volume of the other. The fractions are normalized internally. In contrast, `extracellular_fraction=0.2` means that 20% of the simulated volume is extracellular, with the remaining 80% intracellular.

## Command-line usage

The CLI can create a liminal geometry file using the repeatable `--geometry` option:

```bash
mcmr geometry create liminal cells.txt \
    --extracellular-fraction 0.2 \
    --geometry 0.7 sphere.json \
    --geometry 0.3 cell.swc
```

This creates a text file "cells.txt" with:

```text
liminal 0.2
0.7 sphere.json
0.3 cell.swc
```

The first line identifies the file as a liminal geometry and sets the extracellular volume fraction. Each following line contains a compartment-specific number fraction and a child geometry filename. Child paths are interpreted relative to the liminal file. Child files may be JSON, PLY, or SWC files, and their formats are detected from their contents.

The resulting file can be used directly in a simulation:

```bash
mcmr run cells.txt sequence.seq -N 10000 -o signal.csv
```

The `--geometry` option can be supplied multiple times. Its fractions describe the relative numbers of the listed cell types; they should not be confused with the extracellular volume fraction on the first line.
