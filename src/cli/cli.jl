"""
Defines the command line interface to MCMRSimulator.jl (`mcmr`).

Functions:
- [`run_main`](@ref)
"""
module CLI
include("run.jl")
include("geometry.jl")
include("sequence.jl")

import ArgParse
import Markdown
import .Run
import .Geometry
import .Sequence

function run_main(args=ARGS; kwargs...)
    if length(args) == 0
        println(stderr, "No mcmr command given.\n")
    else
        if args[1] == "run"
            return Run.run_main(args[2:end]; kwargs...)
        elseif args[1] == "geometry"
            return Geometry.run_main(args[2:end]; kwargs...)
        elseif args[1] == "sequence"
            return Sequence.run_main(args[2:end]; kwargs...)
        else
            println(stderr, "Invalid mcmr command $(args[1]) given.\n")
        end
    end
    println(stderr, "usage: mcmr {run/geometry/sequence}")
    return Cint(1)
end

function cmdline_debug_handler(settings::ArgParse.ArgParseSettings, err, err_code::Int = 1)
    println(stderr, err.text)
    println(stderr, ArgParse.usage_string(settings))
    error()
end

"""
    run_main_test(cmd::AbstractString)

Returns the stdout and stderr produced by running [`run_main`](@ref) on the given command.
This is used for testing of the command line interface
"""
function run_main_test(cmd::AbstractString)
    original_stderr = stderr
    err_rd, _ = redirect_stderr()
    err_text = @async read(err_rd, String)
    original_stdout = stdout
    out_rd, _ = redirect_stdout()
    out_text = @async read(out_rd, String)
    prev_interactive = Base.is_interactive
    Base.is_interactive = false

    try
        try
            run_main(split(cmd), exc_handler=cmdline_debug_handler, exit_after_help=false)
        catch err
            if hasproperty(err, :msg) && length(err.msg) > 0
                rethrow(err)
            end
        end
    finally
        close(err_rd)
        redirect_stderr(original_stderr)
        close(out_rd)
        redirect_stdout(original_stdout)
        Base.is_interactive = prev_interactive
    end

    return (fetch(out_text), fetch(err_text))
end

"""
    run_main_test(cmd::AbstractString)

Returns markdown with the stdout and stderr produced by running [`run_main`](@ref) on the given command.
This is used for testing of the command line interface
"""
function run_main_docs(cmd::AbstractString)
    (out, err) = run_main_test(cmd)
    text = ""
    if length(out) > 0
        text = text * "```stdout\n$(out)\n```\n"
    end
    if length(err) > 0
        text = text * "```stderr\n$(err)\n```\n"
    end
    return Markdown.parse(text)
end

end