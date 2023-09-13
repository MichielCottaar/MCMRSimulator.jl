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
        println("No mcmr command given.\n")
    else
        if args[1] == "run"
            return Run.run_main(args[2:end])
        elseif args[1] == "geometry"
            return Geometry.run_main(args[2:end])
        else
            println("Invalid mcmr command $(args[1]) given.\n")
        end
    end
    println("usage: mcmr {run/geometry}")
    return Cint(1)
end
end