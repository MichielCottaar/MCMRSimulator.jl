"""
Defines the command line interface to `mcmr run`.
"""
module Run

import ArgParse: ArgParseSettings, @add_arg_table!, add_arg_group!, parse_args
import DataFrames: DataFrame
import CSV
import ...Geometries.User.JSON: read_geometry
import ...Sequences.PulseQ: read_pulseq
import ...Simulations: Simulation
import ...Spins: Snapshot, BoundingBox, longitudinal, transverse, phase, orientation, position
import ...Evolve: readout


"""Add simulation definition parameters to an argument parser."""
function add_simulation_definition!(parser)
    add_arg_group!(parser, "Define the simulation parameters", :simulation)
    @add_arg_table! parser begin
        "geometry"
            required = true
            help = "JSON file describing the spatial configuration of any obstructions as well as biophysical properties associated with those obstructions. Can be generated using `mcmr geometry`. Alternatively, a mesh file can be provided."
        "sequence"
            nargs = '*'
            help = "One of more pulseq .seq files describing the sequences to be run."
        "--R1"
            help = "Longitudinal relaxation in 1/ms. This relaxation rate will at the very least be applied to free, extra-cellular spins. It might be overriden in the 'geometry' for bound spins or spins inside any obstructions."
            arg_type = Float64
            default = 0.
        "--R2"
            help = "Transverse relaxation in 1/ms. This relaxation rate will at the very least be applied to free, extra-cellular spins. It might be overriden in the 'geometry' for bound spins or spins inside any obstructions."
            arg_type = Float64
            default = 0.
    end
end

"""Add output flags to an argument parser."""
function add_output_flags!(parser)
    add_arg_group!(parser, "Output flags. At least one is required", :output, required=true)
    @add_arg_table! parser begin
        "-o" "--output-signal"
            help = "Writes the total signal at the readouts to this file as a comma-separated value (CSV) table."
        "--output-snapshot"
            help = "Writes the state of all the spins at the readouts to this file as a comma-separated value (CSV) table."
    end
end

"""Add readout flags to an argument parser."""
function add_readout_flags!(parser)
    add_arg_group!(parser, "Readout flags. These control when the signal/spin states will be read out", :readout)
    @add_arg_table! parser begin
        "--nTR"
            help = "Acquire the signal provided at the sequence readouts for this many repetition times (TRs). Output will be stored as a CSV file."
            arg_type = Int
            default = 1
        "-T" "--times"
            help = "Acquire the signal at the given times within each TR (in ms). Multiple values can be provided (e.g., '-T 0 10 15.3'). By default, the Readout markers in the sequence will be used instead."
            arg_type = Float64
            nargs = '+'
        "--skip-TR"
            help = "The number of repetition times the simulation will run before starting to acquire data."
            arg_type = Int
            default = 0
        "--subset"
            help = "Can be provided multiple times. For each time it is provided, the signal will be computed at each readout for a specific subset of spins. Each flag should have 4 values: all/bound/free all/inside/outside <geometry_index> <obstruction_index>."
            nargs = 4
            action = :append_arg
    end
end

"""Add initialisation flags to an argument parser."""
function add_init_flags!(parser)
    add_arg_group!(parser, "Initialisation flags. These control the spins initial state", :init)
    @add_arg_table! parser begin
        "--N"
            help = "Number of spins to simulate. Ignored if --init is set."
            arg_type = Int
            default = 10000
        "--voxel-size"
            help = "Size of the voxel (in mm) over which the initial spins are spread."
            arg_type = Float64
            default = 1.
        "--longitudinal"
            help = "Initial value of the longitudinal magnetisation for each spin. Note the the equilibrium longitudinal magnetisation for each spin is 1."
            arg_type = Float64
            default = 1.
        "--transverse"
            help = "Initial value of the magnitude of the transverse magnetisation for each spin."
            arg_type = Float64
            default = 0.
        "--phase"
            help = "Initial value of the phase of the transverse magnetisation for each spin in degrees."
            arg_type = Float64
            default = 0.
        "--init"
            help = "Continues the simulation with the spin states stored in a JSON file (produced using the 'mcmr pre-run' command). If used all other initialisation flags (except for --seed) are ignored."
        "--seed"
            help = "Initialisation for random number seed. Supply this to get reproducible results. If --init is also set, this flag will override the seed stored in this initialisation file."
            arg_type = Int
    end
end

"""
    get_parser()

Returns the parser of arguments for `mcmr run`
"""
function get_parser()
    parser = ArgParseSettings(prog="mcmr run", description="Runs a Monte Carlo simulation of the MRI signal evolution for spins interacting with the geometry.")
    add_simulation_definition!(parser)
    add_output_flags!(parser)
    add_readout_flags!(parser)
    add_init_flags!(parser)
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
    if !isnothing(args["init"])
        error("Reading snapshots not yet implemented!")
    else
        bb = BoundingBox{3}(args["voxel-size"]/2)
        snap = Snapshot(args["N"], simulation, bb; longitudinal=args["longitudinal"], transverse=args["transverse"])
    end
    as_snapshot = !isnothing(args["output_snapshot"])
    result = readout(snap, simulation, args["times"]; skip_TR=args["skip_TR"], nTR=args["nTR"], noflatten=true, return_snapshot=as_snapshot)

    # convert to tabular format
    if !isnothing(args["output_signal"])
        df_list = []
        for index in eachindex(:IndexCartesian, result)
            if as_snapshot
                value = SpinOrientation(result[index])
            else
                value = result[index]
            end
            orient_as_vec = orientation(value)
            push!(df_list, (
                sequence=index[1],
                TR=index[3],
                readout=index[2],
                subset=1,
                nspins=length(snapshot),
                longitudinal=longitudinal(value),
                transverse=transverse(value),
                phase=phase(transverse),
                Sx=orient_as_vec[1],
                Sy=orient_as_vec[2],
            ))
        end
        df = DataFrame(df_list)
        CSV.write(args["output_signal"], df)
    end

    if !isnothing(args["output_snapshot"])
        df_list = []
        for index in eachindex(:IndexCartesian, result)
            snapshot = result[index]
            for (ispin, spin) in enumerate(snapshot)
                orient_as_vec = orientation(spin)
                pos =  position(spin)
                push!(df_list, (
                    sequence=index[1],
                    TR=index[3],
                    readout=index[2],
                    spin=ispin,
                    x=pos[1],
                    y=pos[2],
                    z=pos[3],
                    longitudinal=longitudinal(spin),
                    transverse=transverse(spin),
                    phase=phase(transverse),
                    Sx=orient_as_vec[1],
                    Sy=orient_as_vec[2],
                ))
            end
        end
        df = DataFrame(df_list)
        CSV.write(args["output_snapshot"], df)
    end

    return Cint(0)
end


end