# Agent Instructions

## Project

- This is a single Julia package with `test` and `docs` subprojects in the root Julia 1.12 workspace; use the root `Project.toml` for package work and the subproject environments for tests/docs.
- The package entrypoint is `src/MCMRSimulator.jl`; it wires together simulation, geometry, Pulseq, plotting, and CLI modules. Geometry code is split between user-facing definitions under `src/geometries/user/` and internal collision/field implementations under `src/geometries/internal/`.
- The CLI entrypoint is `MCMRSimulator.CLI.@main`, with `run` and `geometry` subcommands. The `mcmr` app is declared in `Project.toml` and installed through `Pkg.Apps.add`/`pkg> app add`.

## Commands

- Instantiate dependencies with `julia --project -e 'using Pkg; Pkg.instantiate()'`.
- Run all tests with `julia --project -e 'using Pkg; Pkg.test()'`.
- Test-only dependencies are declared in `test/Project.toml`; the workspace resolves them into the ignored root `Manifest.toml`.
- Run selected test files through `Pkg.test` using `test_args`, for example `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["collisions", "evolve"])'`; use `julia --project -e 'using Pkg; Pkg.test("MCMRSimulator", test_args=["no-plots"])'` to skip the visual `plots` suite.
- CI installs `xvfb` and adds MRIBuilder from its GitLab repository before testing; use `xvfb-run julia --project=@. -e 'using Pkg; Pkg.add(url="https://git.fmrib.ox.ac.uk/ndcn0236/mribuilder.jl.git"); Pkg.test(coverage=true)'` when reproducing CI rather than assuming a plain test run has identical setup.
- Build documentation with `julia --project=docs -e 'using Pkg; Pkg.instantiate(); include("docs/make.jl")'`. `docs/Project.toml` is a workspace subproject, and `docs/make.jl` temporarily generates `docs/src/installation.md` and removes it afterward; do not edit that generated file.

## Change Boundaries

- Update or add tests under `test/` when changing behavior. Visual tests and some CLI/plot paths may need a display or `xvfb`.
- Keep dependency changes synchronized with the appropriate `Project.toml` and `Manifest.toml`; the root `.gitignore` ignores the root manifest even though a local checkout may contain one.
- User-facing MRI acquisition/sequence protocol features belong in the MRIBuilder project, not this package; this repository focuses on simulation, geometry, diffusion, magnetization, and related signal evolution.
