"""
Defines the command line interface to MCMRSimulator.jl (`mcmr`).

Functions:
- [`run_main`](@ref)
"""
module CLI
include("run.jl")
include("geometry.jl")

import .Run
import .Geometry

function run_main(args=ARGS)
    if length(args) == 0
        error("No mcmr command given. Available commands: run, geometry")
    end
    if args[1] == "run"
        Run.run_main(args[2:end])
    elseif args[1] == "geometry"
        Geometry.run_main(args[2:end])
    else
        error("Invalid mcmr command $(args[1]) given. Expected one of: run, geometry")
    end
end
end