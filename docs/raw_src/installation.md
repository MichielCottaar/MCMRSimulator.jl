# [Installation](@id installation)
MCMRSimulator is an application written in the [Julia](https://julialang.org) language.
You can run simulations either directly from the Julia REPL, in a [Jupyter notebook](https://jupyter.org), or using the command line interface.

For now, it is only possible to install MCMRSimulator using the built-in Julia package manager. 
In the future, we will also provide docker and singularity images to make it possible to run the command line interface without installing julia.
# [Global Julia installation](@id global_julia)
1. First install julia from the [official website](https://julialang.org/downloads/).
2. Choose a directory, where you will install MCMRSimulator. If you want to install a single version of MCMRSimulator, this might just be any folder in your home directory. If you want to associate the MCMRSimulator with a specific project, you might want to select a directory within that project. We will refer to this chosen directory below as "<project_dir>".
2. Start the julia REPL in a terminal (`$ julia --project=<project_dir>`).
3. Enter the package manager by pressing "]"
   - Install MCMRSimulator.jl using `pkg> add https://git.fmrib.ox.ac.uk/ndcn0236/mcmrsimulator.jl.git{install_version}`.
   - Optionally, install one of the [Makie backends](https://makie.juliaplots.org/stable/documentation/backends/) for plotting (e.g., `pkg> add CairoMakie`).
   - If you want to use a Jupyter notebook, you will also have to install `IJulia`. You can find instructions to do so at https://github.com/JuliaLang/IJulia.jl.
   - Press "\[backspace\]" to leave the package manager.

## Running MCMRSimulator
After this installation process, you can run MCMRSimulator in one of the following ways:
- *Julia REPL*: Start the REPL in a terminal by typing `$ julia --project=<project_dir>`. Afterwards type `using MCMRSimulator` to import the simulator. You can now follow the steps in the [MCMRSimulator tutorial using Julia](@ref tutorial_julia).
- *Jupyter notebook*: Make sure that you install `IJulia` using the instructions above. This will allow you to start a notebook in jupyter running in Julia. Within this notebook, you can follow the steps in the [MCMRSimulator tutorial using Julia](@ref tutorial_julia).
- *Command line interface*: You can now run the command line interface using `julia --project=<project_dir> -e 'import MCMRSimulator.CLI: run_main; run_main()' -- [args]`. This is a lot to type, so I would recommend adding an alias for this to your ".bashrc" file (or equivalent): `alias mcmr="julia --project=<project_dir> -e 'import MCMRSimulator.CLI: run_main; run_main()' -- "`. With this alias set up, you can now follow [the command line tutorial](@ref tutorial_cli)

## Upgrading MCMRSimulator
1. Start the julia REPL again in a terminal (`$ julia --project=<project_dir>`)
2. Enter the package manager by pressing "]"
3. Update all installed packages using by typing `update` and pressing enter (`pkg> update`).

## Sharing your MCMRSimulator installation
To share the exact environment used by your installation of MCMRSimulator, simply go to the `<project_dir>` directory and locate the files named "Project.toml" and "Manifest.toml". Transfer these files to any other computer, to ensure that they install the exact same version of all Julia packages used (see https://pkgdocs.julialang.org/v1/environments/ for more details).
