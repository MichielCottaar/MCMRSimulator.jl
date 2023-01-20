# MCMRSimulator

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/dev)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/stable)
[![Build Status](https://git.fmrib.ox.ac.uk/ndcn0236/MRSimulator.jl/badges/main/pipeline.svg)](https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/pipelines)
[![Coverage](https://git.fmrib.ox.ac.uk/ndcn0236/MRSimulator.jl/badges/main/coverage.svg)](https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/commits/main)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7318657.svg)](https://doi.org/10.5281/zenodo.7318657)


MRI comes in a wide variety of different modalities, each of which provides its own window into tissue microscture.
Because of these different sensitivities to the tissue, each MRI modality comes with its own microsctural model.
This simulator aims to combine all of these models to produce a single unified tool to model the effect of microstructure on the MRI signal evolution.

## User documentation
The latest documentation can be found [here](https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/dev).

## Citing MCMRSimulator.jl
This software can be cited using the information in [CITATION.cff].

```
@software{cottaar_michiel_2022_7318657,
  author       = {Cottaar, Michiel},
  title        = {MCMRSimulator.jl},
  month        = nov,
  year         = 2022,
  publisher    = {Zenodo},
  version      = {V0.4},
  doi          = {10.5281/zenodo.7318657},
  url          = {https://doi.org/10.5281/zenodo.7318657}
}
```

## Developer documentation
### Release procedure
- Update the version number in "Project.toml" and "README.md" citation section.
- Update the "CHANGELOG.md"
  - Check `[Unreleased]` link for any missing additions to the Changelog
  - Add line with `## [v<version number>]` just below `## [Ureleased]`
  - Add new link at bottom: `[v<version number>]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v<previour version>...v<version_number>`
  - Update unreleased link at bottom with new version number: `[Unreleased]: https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/-/compare/v<previour version>...main`
- Login into [zenodo](https://doi.org/10.5281/zenodo.7318656)
  - In the MCMRSimulator.jl repository click "New version"
  - Click "Reserve doi"
  - Add new citation information to CITATION.cff
    ```
      - description: "This is the archived snapshot of version <version number> of MCMRSimulator.jl"
        type: doi
        value: <reserved doi>
    ```
  - Keep this page open
- Commit changes
- Add tag `v<version number>
- git push
- Create snapshot on gitlab
  - Upload spanshot to zenodo
  - Update version number
- Check if documentation updated correctly