# Extracellular Space Feature Plan

## Goal

Extend SWC and other cell geometries so a simulation can represent both:

- Intracellular spins diffusing inside individual, locally defined cell geometries.
- Extracellular spins diffusing through a simplified, unpositioned liminal space.

Cell geometries do not need globally defined relative positions. A cell is placed in
the laboratory frame only when a spin is assigned to it. Each initially sampled
intracellular spin receives its own independently drawn cell offset, even when
multiple spins use the same cell geometry template. A cell selected by a liminal
encounter also receives an offset unique to that spin.

## User Model

The user provides:

- A collection of cell geometry templates.
- One normalized volume fraction per cell geometry.
- An extracellular volume fraction `f_ext`.
- Optionally, a laboratory-frame bounding box.

The normalized intracellular fractions are scaled as:

```text
f_i = (1 - f_ext) * w_i / sum(w)
```

where `w_i` is the user-provided normalized fraction for cell geometry `i`.

The extracellular fraction is:

```text
f_ext
```

The user-provided bounding box is used for initial cell offsets and liminal spin
positions. If no bounding box is supplied, use the centered 1 mm^3 box currently
represented by `BoundingBox(500)`.

## Physical Model

### Intracellular diffusion

Intracellular diffusion uses the existing geometry collision machinery:

- The spin has an active cell geometry index.
- The spin has an active cell offset.
- Collision detection is performed in the active cell's local coordinates.
- The local collision result is transformed back to laboratory coordinates.
- Existing surface permeability, surface relaxation, sticking, and cached inside
  indices remain applicable.

The active cell offset is not a global geometry placement. It is state carried by the
spin, because cells are templates and every spin may use a different placement.

### Liminal diffusion

Liminal spins diffuse freely in laboratory coordinates. They do not run normal
geometry collision detection until an encounter is sampled.

For a travel direction `u`, calculate the inverse mean free path as:

```text
lambda_inv(u) = sum_i f_i * A_projected_i(u) / V_i
```

where:

- `f_i` is the intracellular volume fraction of cell geometry `i`;
- `A_projected_i(u)` is the total outer surface area projected onto `u`;
- `V_i` is the volume of cell geometry `i`.

For a liminal displacement with length `d`, sample an encounter with probability:

```text
p_encounter = 1 - exp(-lambda_inv(u) * d)
```

If no encounter occurs, the spin moves to the proposed position.

If an encounter occurs:

1. Select a cell geometry and surface sample using projected-area weights.
2. Select a local outer-surface sample with weight proportional to:
   `surface_weight * abs(normal ⋅ u)`.
3. Translate the selected cell so that the sampled local surface position coincides
   with the encounter position.
4. Apply the selected surface's permeability.

An impermeable encounter reflects the spin and returns it to liminal diffusion. A
permeable encounter enters the selected shifted cell and continues with normal
intracellular collision handling.

### Detailed balance

The outer encounter rate must match the intracellular outward flux for each outer
surface element. The equilibrium test must verify that the extracellular fraction
remains stable over long diffusion runs, within statistical tolerance.

## Surface Preprocessing

Each cell geometry needs a precomputed outer-surface sample library.

Use the existing `random_surface_spins` or lower-level
`random_surface_positions` machinery as the starting point. The preprocessing must:

- Sample a large, configurable number of surface points.
- Remove samples that lie inside another part of the same cell geometry.
- Retain only samples on the outer surface.
- Store the local position.
- Store the outward normal.
- Store the surface/property index.
- Store a surface-element weight.

The surface weights must allow the samples to estimate total surface area and
projected surface area. If the existing sampler does not expose sufficient area
weights, add a preprocessing-specific sampler or attach a uniform quadrature weight
based on the represented surface area and number of retained samples.

The preprocessing result should contain at least:

```julia
OuterSurfaceSamples(
    positions,
    normals,
    weights,
    property_indices,
    total_area,
    volume,
)
```

For a direction `u`, projected area is estimated by:

```julia
sum(sample.weight * abs(dot(sample.normal, u)) for sample in samples)
```

The preprocessing must be deterministic when supplied with a seeded RNG, so tests
can compare projected-area estimates and encounter selection behavior.

## Proposed Representation

### User-facing collection

Introduce a user-facing collection type, with a name chosen to match the package's
existing geometry naming conventions. It should contain:

- Cell geometry templates.
- Normalized volume fractions.
- Extracellular volume fraction.
- Optional bounding box.
- Optional surface preprocessing size/seed settings.

Possible conceptual API:

```julia
cells = CellCollection(
    geometries=[read_swc("cell_a.swc"), read_swc("cell_b.swc")],
    volume_fraction=[1.0, 2.0],
    extracellular_fraction=0.2,
    bounding_box=BoundingBox(500),
)
```

The final public name and exact keyword names should follow the existing geometry
API after the internal prototype is working.

### Internal fixed representation

Represent the collection as a vector of independently fixed cell templates. Each
entry should contain:

- A fixed physical geometry.
- Its fixed property metadata.
- Its normalized/normalized-to-total volume fraction.
- Its outer-surface samples.
- Its local volume and projected-area data.

Do not flatten all cell templates into one globally positioned `GeometryTuple`.
Flattening would lose the distinction between cell-local coordinates and per-spin
cell offsets.

The existing property representation can still use tuple-like indexing for each
cell's local properties. Extracellular relaxation and off-resonance should use the
global simulation properties, not cell surface properties.

## Spin State

The current spin state is specialized around `Nothing` or `Reflection`. Extend it
with an explicit movement state rather than encoding liminal behavior as a fake
reflection.

The state needs to distinguish:

```text
LiminalState
    laboratory position
    no active cell

IntracellularState
    active cell index
    active cell offset
    existing reflection/bound state, if any
    cached local inside indices
```

Possible implementation strategies:

1. Add a separate state field to `Spin` while retaining `Reflection` for surface
   collision state.
2. Generalize the existing reflection type parameter into a movement-state type.

Prefer the smallest design that preserves existing `Spin` construction, deepcopy,
snapshot slicing, and multi-sequence behavior. Avoid adding backward compatibility
unless existing persisted or external spin state requires it.

The state must survive:

- `deepcopy`.
- `Snapshot` construction and slicing.
- `Snapshot(snap, nsequences)`.
- `_to_snapshot` and `_constrain_snapshot`.
- `readout(..., return_snapshot=true)`.
- Passing a returned snapshot into another `evolve` call.

## Initial Spin Sampling

When constructing `Snapshot(nspins, simulation, bounding_box)` for the new
collection:

1. Resolve the bounding box.
   - Use the explicit box when supplied.
   - Otherwise use the centered 1 mm^3 box.
2. Allocate extracellular spins according to `f_ext`.
3. Allocate intracellular spins according to the scaled cell fractions `f_i`.
4. Draw liminal spin positions uniformly in the laboratory bounding box.
5. For each cell geometry receiving intracellular spins:
   - sample local positions inside the cell geometry;
   - draw one independent random laboratory offset for each sampled spin;
   - apply that offset to the sampled spin's cell geometry;
   - initialize their active cell index and local inside cache.

The exact integer allocation should be deterministic for a fixed RNG. Use a clear
rounding rule and test that the resulting fraction is within one spin of the target
fraction, or use a binomial allocation if statistical initialization is preferred.

The existing ordinary-geometry snapshot behavior must remain unchanged.

## Evolution State Transitions

Refactor the core movement loop in `draw_step!` into state-specific operations.

### Intracellular step

1. Draw the normal Gaussian displacement in cell-local coordinates.
2. Detect the closest collision against the active cell template.
3. Apply relaxation over the traveled portion.
4. Apply surface relaxation and permeability.
5. On a reflected collision, remain intracellular in the same cell.
6. On a permeable outward crossing, transform to laboratory coordinates and enter
   `LiminalState`.
7. On a permeable inward crossing from an encounter, enter the selected cell state.

### Liminal step

1. Draw a Gaussian laboratory-frame displacement.
2. Compute its direction and length.
3. Compute projected inverse mean free path from all cell templates.
4. Sample an encounter distance/probability.
5. If there is no encounter, finish the step in liminal space.
6. If there is an encounter, place the selected cell surface sample at the encounter
   position.
7. Apply surface permeability.
8. If reflected, continue the remaining displacement in liminal space.
9. If permeable, enter the selected shifted cell and process the remaining
   displacement intracellularly.

The implementation must preserve displacement and time accounting when an encounter
occurs partway through a timestep. The existing reflection bookkeeping can be used
as a reference, but liminal encounters need a separate path because they are sampled
probabilistically rather than found by ray intersection.

## Properties

For an encounter, use the selected cell surface sample's property index to resolve
membrane permeability and any cell-specific membrane behavior.

For liminal movement:

- Use global simulation `R1`/`R2` and off-resonance parameters.
- Do not apply a cell surface relaxation unless an encounter actually occurs.
- Do not create a bound/stuck state for a liminal spin unless the selected membrane
  interaction explicitly supports it.

The timestep controller must account for the new encounter probability. Add a
collection-level maximum inverse mean free path or equivalent timestep constraint if
needed to keep encounter probabilities numerically well resolved.

## Implementation Steps

### Step 1: Document and lock down existing interfaces

- Confirm the existing volume calculation for each fixed cell geometry.
- Confirm how `random_surface_positions` exposes normals, property indices, and
  surface metadata.
- Confirm how outer-surface filtering works for overlapping components.
- Add helper tests before changing `Spin` or `draw_step!`.

### Step 2: Build the fixed cell collection

- Add the user configuration type and validation.
- Validate nonempty geometries, positive normalized fractions, and
  `0 <= extracellular_fraction <= 1`.
- Fix each geometry independently.
- Compute the scaled intracellular fractions.
- Resolve the collection bounding-box behavior.
- Keep this representation separate from ordinary fixed geometry.

### Step 3: Implement outer-surface preprocessing

- Generate the large local surface sample library for each cell.
- Filter internal surface samples.
- Attach normals, property indices, and area weights.
- Add projected-area estimation for arbitrary directions.
- Add tests for spheres, finite cylinders, overlapping SWC components, and mixed cell
  templates.

### Step 4: Add explicit cell/liminal spin state

- Extend `Spin` with the smallest suitable state representation.
- Update constructors and type inference.
- Update deepcopy and snapshot conversion paths.
- Add state round-trip tests without evolution.

### Step 5: Implement initial population sampling

- Allocate intracellular and extracellular spins.
- Sample liminal positions in the laboratory bounding box.
- Sample intracellular local positions.
- Apply one independent random offset per initially sampled intracellular spin.
- Verify initial population fractions and cell membership.

### Step 6: Refactor evolution into state-specific paths

- Isolate the existing intracellular path first without changing behavior.
- Add liminal free-diffusion movement.
- Add projected encounter probability.
- Add weighted cell and surface-sample selection.
- Add per-spin cell offset placement.
- Handle permeable entry and impermeable liminal reflection.
- Preserve remaining timestep displacement across transitions.

### Step 7: Add equilibrium and detailed-balance validation

- Run many spins for a sufficiently long diffusion period.
- Verify intracellular and extracellular fractions remain close to their requested
  values.
- Test a single simple sphere analytically before testing SWC cells.
- Test multiple cell geometries with different volumes and surface areas.
- Test both isotropic and directionally biased displacements.

### Step 8: Integrate user-facing workflows

- Add constructors and exports.
- Add documentation and examples using one or more SWC files.
- Add CLI support only if the feature belongs in this package's existing CLI scope.
- Update `CHANGELOG.md`.
- Build the documentation.

## Test Plan

### Geometry and preprocessing

- Outer surface filtering removes internal surfaces.
- Surface sample normals point outward.
- Projected area of a sphere agrees with its analytic projection.
- Projected area estimates converge as sample count increases.
- Finite-cylinder side and cap samples are handled correctly.
- Mixed SWC geometries preserve per-cell property indices.

### Initialization

- Explicit bounding boxes are honored.
- Missing bounding boxes use the centered 1 mm^3 default.
- Liminal spins are uniformly distributed in the box.
- Initial intracellular/extracellular counts match the requested fractions.
- Intracellular spins using the same geometry template receive independent offsets.
- Different spins receive independent cell placements.

### Transitions

- Intracellular permeable exit produces a liminal spin.
- Liminal motion does not collide with geometry unless an encounter is sampled.
- Encounter selection follows projected-area weighting.
- Encountered cells are shifted so the sampled surface is at the encounter point.
- Impermeable encounters reflect and remain liminal.
- Permeable encounters enter the selected cell.
- Remaining displacement is handled correctly after an encounter.

### Equilibrium

- A single sphere preserves the requested extracellular fraction.
- Opposing and overlapping SWC components do not cause bounce loops.
- Multiple geometries preserve the normalized intracellular fractions.
- Long runs do not systematically drift toward intracellular or extracellular space.

### State persistence

- Deepcopy preserves active cell, offset, liminal state, and inside cache.
- Readout snapshots preserve state.
- Evolving a returned snapshot continues from the same state.
- Multi-sequence conversion preserves movement state.
- Snapshot slicing does not share mutable state unexpectedly.

## Performance Considerations

- Surface samples should be generated once per fixed cell template, not per spin.
- Projected-area weights should be cached or computed using allocation-free loops.
- Encounter sampling should avoid scanning millions of samples for every spin if a
  direction-dependent alias table or hierarchical sampler can be introduced later.
- Start with a straightforward weighted scan for correctness, then benchmark.
- Make the preprocessing sample count configurable for production versus tests.

## Risks and Safeguards

- **Surface sampling bias:** retain explicit area weights and test analytic shapes.
- **Coincident surfaces:** filter internal surfaces and preserve the existing
  coincident-cap safeguards for SWC cylinders.
- **State corruption:** centralize transitions and add snapshot round-trip tests.
- **Incorrect detailed balance:** test population fractions over long runs, not only
  individual transitions.
- **Arbitrary gradient phase:** all initial liminal and cell offsets must be sampled
  in the same laboratory bounding box so gradient-dependent signal has a defined
  coordinate frame.
- **Large preprocessing cost:** expose sample count and deterministic RNG settings.

## Acceptance Criteria

The feature is ready when:

- Individual SWC cells can be used as local intracellular templates.
- Initial spins include both intracellular and liminal populations.
- Liminal encounters select cells and surface points according to projected area.
- Permeable and impermeable encounters produce the correct state transitions.
- The requested extracellular fraction remains stable in long diffusion tests.
- Existing non-extracellular tests remain unchanged and pass.
- Documentation contains a reproducible example using SWC cell geometries.
