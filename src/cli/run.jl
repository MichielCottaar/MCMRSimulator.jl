"""
Defines the command line interface to `mcmr run`.
"""
module Run

import ArgParse: ArgParseSettings, @add_arg_table!, add_arg_group!, parse_args
import ...Geometries.User.JSON: read_geometry
import ...Sequences.PulseQ: read_pulseq
import ...Simulations: Simulation
import ...Spins: Snapshot, BoundingBox

"""
    get_parser()

Returns the parser of arguments for `mcmr run`
"""
function get_parser()
    parser = ArgParseSettings(prog="mcmr run", description="Runs a Monte Carlo simulation of the MRI signal evolution for spins interacting with the geometry.")
    add_arg_group!(parser, "output options. Exactly one of these options should be supplied.", :output, false, required=true, exclusive=true)
    @add_arg_table! parser begin
        "geometry"
            required = true
            help = "JSON file describing the spatial configuration of any obstructions as well as biophysical proeprties associated with those obstructions. Can be generated using `mcmr geometry`. Alternatively, a mesh file can be provided."
        "sequence"
            required = true
            nargs = '+'
            help = "One of more pulseq .seq files describing the sequences to be run."
        "-o" "--output"
            help = "The output filename. The format of the output file depends on which of the other output options are provided."
            required = true
        "--nTR"
            help = "Acquire the signal provided at the sequence readouts for this many repetition times (TRs). Output will be stored as a CSV file."
            arg_type = Int
            group = :output
        "--TE"
            help = "Acquire the signal at the given times (in ms). Output will be stored as a CSV file."
            arg_type = Float64
            nargs = '+'
            group = :output
        "--snapshot"
            help = "Output the full state of the simulation at given time (in ms). The output will be stored in a JSON file."
            arg_type = Float64
            group = :output
        "-N"
            help = "Number of spins to simulate (default 10,000)."
            arg_type = Int
            default = 10000
        "--voxel-size"
            help = "Size of the voxel (in mm) over which the initial spins are spread (default: 1 mm)."
            arg_type = Float64
            default = 1.
        "--seed"
            help = "Initialisation for random number seed. Supply this to get reproducible results."
            arg_type = Int
        "--init"
            help = "Continues the simulation with the spin states stored in a JSON file (produced using the `--snapshot` flag)."
        "--R1"
            help = "Longitudinal relaxation rate (i.e., 1/T1 in 1/ms). This might be overriden for bound spins or within certain obstructions by `geometry`."
            default = 0.
            arg_type = Float64
        "--R2"
            help = "Transverse relaxation rate (i.e., 1/T2 in 1/ms). This might be overriden for bound spins or within certain obstructions by `geometry`."
            default = 0.
            arg_type = Float64
        "--skip-TR"
            help = "The number of repetition times the simulation will run before starting to acquire data. The time provided in the output options (`--TE` or `--snapshot`) are relative to the end of this equilibriation period."
            arg_type = Int
            default=0
        "--subset"
            help = "Can be provided multiple times. For each time it is provided, the signal will be added for a specific subset of spins. Each flag should have 4 values: all/bound/free all/inside/outside <geometry_index> <obstruction_index>."
            nargs = 4
            action = :append_arg
            
    end
    return parser
end

"""
    run_main([arguments])

Runs the `mcmr run` command line interface.
The supplied arguments can be provided as a sequence of strings (as provided by the terminal)
or as a dictionary (as provided by `ArgParse` after parsing).
By default it is set to `ARGS`.
"""
function run_main(args=ARGS::AbstractVector[<:AbstractString])
    parser = get_parser()
    run_main(parse_args(args, parser))
end

run_main(::Nothing) = Cint(1)

function run_main(args::Dict{<:AbstractString, <:Any})
    geometry = read_geometry(args["geometry"])
    sequences = read_pulseq.(args["sequence"])

    simulation = Simulation(sequences; geometry=geometry, R1=args["R1"], R2=args["R2"])

    # initial state
    if !isnothing(args["snapshot"])
        error("Reading snapshots not yet implemented!")
    else
        bb = BoundingBox(args["voxel-size"]/2)
        snap = Snapshot(args["N"], simulation, bb)
    end
    return Cint(0)
end


end