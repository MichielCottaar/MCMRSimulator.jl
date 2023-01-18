# Changelog
All notable changes to MCMRSimulator.jl will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).
## [Unreleased]
### Added
- Support for arbitrary gradient durations and diffusion times in `dwi` thanks to Zhiyu
- Flexible timesteps that dynamically adjust to the sequence(s)
### Fixed
- Magnetisation transfer rate is now corrected for the dependence of the collision rate on the timestep
- Running `signal` for a `Simulation` with multiple sequences used to crash. This has now been fixed
- Fix bouncing between repeating walls (used to get escapees)
### Documentation
- Expanded README intro
- Added citation instructions


[Unreleased]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v0.3.0...main