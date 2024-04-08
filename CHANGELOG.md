# Changelog
All notable changes to MCMRSimulator.jl will be documented here.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
## [Unreleased]
### Added
- `mcmr sequence plot` can now be used to plot sequence diagrams from the command line.
### Changed
- `random_positions_radii` now has an additional step of repulsion between the cylinders/spheres. This ensures a much smoother density distribution in the output.
- The plotting in MCMRSimulator has been refactored:
    - `plot_sequence` produces an actual sequence diagram
    - More keywords can be passed on to `plot_geometry`, `plot_off_resonance`, `plot_snapshot`, `plot_sequence`, and `plot_trajectory`. The documentation of these functions now reflects these new keywords.
    - The plotting library `Makie` is now an optional dependency. This greatly speeds up loading time, if plotting is not required.
    - `plot_geometry3d`, `plot_trajectory3d`, and `plot_snapshot3d` have been removed. Used `plot_geometry`, `plot_trajectory`, and `plot_snapshot` (or just `plot`) instead.
    - Plotting of individual snapshots can now be done using `Makie.plot([plot_plane, ]snapshot)`. You can decide whether to plot the spins as points (default), dyads (set `kind=:dyad`), or the total magnetisation as an image (set `kind=:image`). The length of the dyads is set using the `lengthscale` keyword (as in `Makie.arrows`), instead of the `dyadlength` keyword (as previous).
## [v0.9.0]
### Added
- Text files can be used to set geometry positions/radii/etc from the command line.
- New obstruction type: `BendyCylinder`. This cylinder can follow an arbitrary path (defined using splines) and vary in diameter.
- 3D plotting for meshes (and `BendyCylinder`) is now available. Just call `plot(mesh)` or `plot_geometry3d(mesh)`.
### Changed
- Multiple values passed on the `mcmr geometry` on the command line should now be separated by commas and semi-colons rather than spaces
- The `save_memory` option has been removed from meshes, because it did not work properly.
### Fixed
- Fixed occasional leakage through `Walls` that are exactly at the edge of a grid cell or repeat.
- Creation of large meshes has been sped up (and no longer leads to `StackOverflow` errors for very large meshes).
- Off-resonance fields generated from meshes are now correctly modelled using magnetic dipoles rather than magnetic monopoles.
- Size calculation of meshes has been fixed (used in the calculation of the maximum timestep).
- Determining the inside of meshes has been fixed.
- Fixed printing of arrays of `Snapshot`.
### Changed
- CSV output files now contain the name of the sequence file (CLI)
## [v0.8.0]
### Added
- New `mcmr run --bvecs` flag and `rotate_bvecs` function that apply a rotation to all diffusion-weighted gradeints in an MRI sequences
- New `mcmr sequence <name> --scanner` flag to select one of the predefined scanners.
### Changed
- Crusher gradients have been added to the pre-defined sequences (`gradient_echo`, `spin_echo`, and `dwi`).
- For a 3D geometry (i.e., mesh or spheres) if the orientation is set as a vector, this vector represents the new orientation of the geometry x-axis (previously it was the z-axis).
### Fixed
- Greatly speed up generation of Snapshots for simulations with large voxels.
- Playing out the sequence occassionally broke down after multiple repetition times. This has been fixed.
- The `--seed` flag is now actually used in `mcmr run` and also available in `mcmr create-random`
- Fixed an issue that caused `mcmr create-random` to crash
## [v0.7.1]
### Fixed
- Fix parsing of `position` flag in CLI.
- Remove `--init` flag from `mcmr run`, because it has not been implemented yet
## [v0.7.0]
### Added
- New command line interface. See the tutorial in the documentation on how to use it.
- `isinside` now works for meshes. This can be used to filter intra- or extra-cellular water as well as to set a different T1/T2/off-resonance in the intra-cellular space.
- Off-resonance fields can now be calculated for meshes. Set `myelin=true` for the mesh to enable this.
- `split_mesh` function that splits a mesh into connected components (and fixes the mesh normals)
### Changed
- Default diffusivity changed from 0 to 3 um^2/ms
- Replaced `readout`, `trajectory`, and `signal` with a single, improved `readout` function. This new function can:
    - return the total signal or snapshots using on `return_snapshot` keyword
    - either return the signal at the `Readout` objects in the sequence or at any times set when calling the `readout` function
    - allow the signal to equilibriate for multiple repetition times before returning an output (using the new `skip_TR` keyword)
    - return the signal for multiple repetition times (using the new `nTR` keyword)
    - This does not affect the `evolve` function, which has not been changed.
- It is made easier to select a specific subset of spins using `get_subset` or by providing one or more `Subset` objects to the `readout` function.
- Geometries are now mutable.
- The `R1`, `R2`, and `off-resonance` values set by the geometry (inside or stuck) are now added to the global values, rather than replacing them.
- Geometries are now generated by calling their plural type (e.g., `Spheres`, `Annuli`, `Mesh`). Geometries can no longer be generated individually using the singular type (e.g., `Sphere`, `Annulus`) or using a lower-case constructor (e.g., `spheres`, `annuli`). The interface of the plural type closely matches that of the constructors, except for:
    - The "MT_fraction" keyword has been renamed "surface_relaxivity".
    - The "positions" keyword in the geometry constructors (e.g., `Walls` or `Annuli`) is now renamed to "position".
    - The radii of `Cylinders` and `Spheres` are now set using the "radius" keyword rather than as positional arguments.
    - The inner and outer radii of `Annuli` are now set using the "inner" and "outer" keywords rather than as positional arguments.
### Fixed
- Fixed the calculation of the off-resonance fields of annuli and cylinders for oblique magnetic field orientations.
- `random_positions_radii` now accepts a variance of zero.
- Sequence plots now work even for zero-amplitude instant gradients or RF pulses.
### Refactor
- The code is now split up into many Julia modules with each file corresponding to a module. This should not affect the package user.
- The geometry module has been rewritten to separate the user interface for setting/updating the geometry from the internal representation of that geometry. This change meant the internal representation could be optimised for running speed. There are some minor changes in the rewrite (see "Changed" above).
- The collision properties (`permeability`, `surface_density`, `dwell_time`, and `surface_relaxivity`) are no longer stored globally in the simulation. Instead, they are integrated into the fixed geometry object.
## [v0.6.0]
### Added
- Spins getting stuck on the surface. The rate of spins getting stuck is controlled by setting the surface density of stuck spins (relative to volume density) and the dwell time of those spins. This could represent surface tension or magnetisation transfer. Stuck spins can be given different MRI relaxation properties. During `Snapshot` initialisation the correct fraction of spins will be initialised as being stuck on the surface if a `Simulation` object is provided. If a `Snapshot` is generated by giving the number of required spins to `readout`, `evolve`, `trajectory`, or `signal`, stuck spins will always be initialised.
- A new sequence builder paradigm using sequence `BuildingBlock` objects as explained in the "Sequence" section in the tutorial.
- New `gradient_echo` and `spin_echo` functions
- v1.4.0 of the `pulseq` sequence is now supported
### Changed
- The spelling of `read_pulseseq` has been corrected to `read_pulseq`.
- `Sequence` constructor now expected all sequence components (pulses, gradients, and readouts) to be passed on as a vector to a single `components` flag rather than separate `pulses` and `gradients` flags.
### Fixed
- A `DomainError` is now thrown when obstructions are cut off by the repeats.
- The example `Scanner` objects were created with the wrong units. Their maximum gradient strengths and slew rates have now been fixed.
- Fixed gradient strengths and RF pulse phases read from `pulseq` sequence files
## [v0.5.0]
### Added
- Finite RF pulses (`RFPulse`) to more realistically model the effect of these pulses.
- Hierchical properties: MRI relaxation parameters and collision parameters (e.g., permeability) can now set at both the global level and for each obstruction.
- Sequences can now be read from the pulseq format (http://pulseq.github.io/) using `read_pulseq`.
- `run_tests.sh` function which will run the tests and produce a coverage that can be visualised in VS Code by the Coverage Gutters plugin.
### Changed
- Timestep calculation has been simplified to only three parts (see `TimeController` or `propose_times` for details), namely (1) Geometry size scale, (2) Gradient size (whether internal due to myelin/iron off-resonance fields or external), and (3) RF pulse maximum rotation. These can be controlled using:
    - `maximum_timestep`: overrides the maximum timestep set by the geometry size scale and the internal gradient size
    - `gradient_size`: sets the error allowed when evaluating the phase evolution due to internal or external gradients (default: 1 degree). This sets an additional constraint on top of the `maximum_timestep` while there are strong gradients.
    - `rf_rotation`: maximum rotation due to RF pulses that can occur in a single timestep (defualt: 1 degree). This will shorten the timestep while RF pulses are active.
- Made units consistent with angles (i.e., phases and flip angles) in degrees, off-resonance fields in kHz, and gradients in kHz/um.
- `random_positions_radii` now allows a minimum and maximum radius to be set. The minimum radius is set to 0.1 by default.
- `time` function has been renamed `get_time`, so as not to conflict with `Base.time`
- Refactored the `Sequence` interface:
    - `RFPulse` -> `InstantRFPulse` (`RFPulse` now refers to finite RF pulses discussed above)
    - `SequenceComponent` -> `InstantComponent` (`Readout` is no longer a sub-type of this abstract type)
- Relaxation parameters can now be set as relaxation times (T1/T2) in addition to as relaxation rates (R1/R2)
- You will no longer have to iterate over sequences in the output of `signal` and `readout` if the simulation is created with a singular sequence object (e.g., `Simulation(sequence, ...)`). Note that this will not affect simulations with a vector of sequence objects, even if the length of that vector is one (e.g., `Simulation([sequence])`).
- Pretty printing for simulations, sequences, spins, and snapshots is much improved.
- Spin magnetisations are now updated between collisions during the movement rather than only between timesteps.

### Fixed
- In the `readout` output, snapshots were included multiple times for the same timepoint if the same readout time was present across the sequences. This has now been fixed, so that for each sequence there will only be a single snapshot per readout time.
- The `nsequenes` flag in the `Spin` constructor now works.
## [v0.4.0]
### Added
- Support for arbitrary gradient durations and diffusion times in `dwi` thanks to Zhiyu
- Flexible timesteps that dynamically adjust to the sequence(s)
- Obstacles can now be made permeable using the `permeablity` keyword.
### Fixed
- Magnetisation transfer rate is now corrected for the dependence of the collision rate on the timestep
- Running `signal` for a `Simulation` with multiple sequences used to crash. This has now been fixed
- Fix bouncing between repeating walls (used to get escapees)
### Documentation
- Expanded README intro
- Added citation instructions


[Unreleased]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.9.0...main
[v0.9.0]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.8.0...v0.9.0
[v0.8.0]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.7.1...v0.8.0
[v0.7.1]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.7.0...v0.7.1
[v0.7.0]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.6.0...v0.7.0
[v0.6.0]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.5.0...v0.6.0
[v0.5.0]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.4.0...v0.5.0
[v0.4.0]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.3.0...v0.4.0