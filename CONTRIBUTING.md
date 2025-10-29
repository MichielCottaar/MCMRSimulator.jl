# Contributing to MCMRSimulator

Thanks for your interest in helping out! The project is currently maintained by Michiel Cottaar, and we are happy to collaborate on improvements that make the simulator more useful to the community.

## Reporting issues

- Please search the [existing issues](https://git.fmrib.ox.ac.uk/ndcn0236/mcmrsimulator.jl/-/issues) before opening a new one.
- When filing a bug report, include:
  - The version of `MCMRSimulator.jl` you are using (`Project.toml` and `Manifest.toml` can help here).
  - A minimal example (code snippets or data) that reproduces the problem.
  - Any relevant logs or stack traces.
- Feature requests can include novel geometries or biophyisical features. MRI acquisition or sequence protocols can be contributed to [MRIBuilder](https://git.fmrib.ox.ac.uk/ndcn0236/mribuilder.jl) instead.
   - *Novel geometries* are shapes that represent a novel implementation of spin collision algorithm or the calculation of the off-resonance fields. New packing algorithms for existing geometries (spheres/cylinders/mesh) are easier to implement as separate packages rather than included in the simulator.
   - *Biophysical features* can include any changes to how the spins move around or how the spin magnetisation is updated.

## Getting started

1. Clone the repository and activate the Julia project:
   ```bash
   git clone https://git.fmrib.ox.ac.uk/ndcn0236/mcmrsimulator.jl
   cd mcmrsimulator.jl
   julia --project -e 'using Pkg; Pkg.instantiate()'
   ```
2. Run the test suite to check everything works locally:
   ```bash
   julia --project -e 'using Pkg; Pkg.test()'
   ```

## Making changes
Please contact Michiel Cottaar before making any substantial changes to the code base, so we can discuss the best way to implement a particular feature and reduce the risk of duplicate efforts.

- Use feature branches for your work and keep commits focused.
- Follow the existing coding style; if in doubt, mimic the surrounding code.
- Update or add tests when fixing bugs or adding functionality.
- Update documentation (including tutorials or docstrings) when user-facing behaviour changes.

## Release procedure
- Check out the branch "`v<major>.<minor>`" (might already exist if this is a patch update)
  - Rebase the main branch into the version branch
- Update the version number in "Project.toml" and "README.md" citation section.
- Update the "CHANGELOG.md"
  - Check `[Unreleased]` link for any missing additions to the Changelog
  - Add line with `## [v<version number>]` just below `## [Unreleased]`
  - Add new link at bottom: `[v<version number>]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v<previous version>...v<version_number>`
  - Update unreleased link at bottom with new version number: `[Unreleased]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v<version number>...main`
- Login into [zenodo](https://doi.org/10.5281/zenodo.7318656)
  - In the MCMRSimulator.jl repository click "New version"
  - Click "Reserve doi" (should already be clicked)
  - Add new citation information to CITATION.cff (do not update doi in README)
    ```
      - description: "This is the archived snapshot of version <version number> of MCMRSimulator.jl"
        type: doi
        value: <reserved doi>
    ```
  - Keep the zenodo page open
- Commit changes
- Add tag `v<version number>`
- git push
- Create and auto-merge merge request from version branch back into main
  -  Ensure "delete branch after merge" is not selected!
- Create release on gitlab (https://git.fmrib.ox.ac.uk/ndcn0236/mcmrsimulator.jl/-/releases)
  - title: `v<version number>`; Release notes: "See CHANGELOG.md for list of changes.";
  - Upload spanshot (.zip) to zenodo
  - Update version number in zenodo
  - Press "Save" and then "Process" buttons at the bottom of the zenodo page.
- After pipelines finish:
  - Check if documentation updated correctly
  - Run `./local_docker_build.sh v<version number>` to build docker ARM64 image
  - Build singularity image on https://cloud.sylabs.io/builder (copy build instructions from last time and update version number in "%labels" section and "From:" statement)