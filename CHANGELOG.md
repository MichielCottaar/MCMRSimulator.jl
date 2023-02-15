# Changelog
All notable changes to MCMRSimulator.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
## [Unreleased]
### Fixed
- A `DomainError` is now thrown when obstructions are cut off by the repeats.
## [0.5.0]
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


[Unreleased]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.5.0...main
[v0.5.0]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.4.0...v0.5.0
[v0.4.0]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.3.0...v0.4.0