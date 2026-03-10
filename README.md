# MCMRSimulator

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/dev)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/stable)
[![Build Status](https://git.fmrib.ox.ac.uk/ndcn0236/MRSimulator.jl/badges/main/pipeline.svg)](https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/pipelines)
[![Coverage](https://git.fmrib.ox.ac.uk/ndcn0236/MRSimulator.jl/badges/main/coverage.svg)](https://git.fmrib.ox.ac.uk/ndcn0236/MCMRSimulator.jl/commits/main)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7318656.svg)](https://doi.org/10.5281/zenodo.7318656)


MRI comes in a wide variety of different modalities, each of which provides its own window into tissue microscture.
Because of these different sensitivities to the tissue, each MRI modality comes with its own microsctural model.
This simulator aims to combine all of these models to produce a single unified tool to model the effect of microstructure on the MRI signal evolution.

## User documentation
The latest documentation can be found [here](https://open.win.ox.ac.uk/pages/ndcn0236/mcmrsimulator.jl/dev).

## Contributing to the simulator
MCMRSimulator is currently maintained by Michiel Cottaar. Bug reports and feature requests are welcome through the [GitLab issue tracker](https://git.fmrib.ox.ac.uk/ndcn0236/mcmrsimulator.jl/-/issues); please check for existing issues before opening a new one. If you would like to collaborate on new features, feel free to start a discussion there so we can work on an implementation together. For more detail, see the [contributor guide](CONTRIBUTING.md).

## Citing MCMRSimulator.jl
If you use the simulator, please cite the following journal article:
```
@article{cottaarMultimodalMonteCarlo2026,
  title = {Multi-Modal {{Monte Carlo MRI}} Simulator of Tissue Microstructure},
  author = {Cottaar, Michiel and Zheng, Zhiyu and Miller, Karla and Tendler, Benjamin C. and Jbabdi, Saad},
  year = 2026,
  month = mar,
  journal = {Imaging Neuroscience},
  issn = {2837-6056},
  doi = {10.1162/IMAG.a.1177},
}
```