"""
Defines command line interface for `mcmr geometry`
"""
module Geometry

import ArgParse: ArgParseSettings, @add_arg_table!, add_arg_table!, parse_args, ArgParseError, usage_string
import StaticArrays: SizedVector
import ...Geometries.User.Obstructions: walls, Walls, spheres, Spheres, cylinders, Cylinders, annuli, Annuli, fields, field_to_docs, Field, ObstructionGroup
import ...Geometries.User.JSON: write_geometry

field_type(::Field{T}) where {T} = T

"""
    get_parser()

Returns the parser of arguments for `mcmr geometry`
"""
function get_parser()
    parser = ArgParseSettings(
        prog="mcmr geometry", 
        description="Various commands to manipulate the geometry used in MCMR simulations",
    )
    @add_arg_table! parser begin
        "create"
            help="Creates a geometry JSON file containing a single set of obstructions (e.g., cylinders, spheres, walls)."
            action=:command
        "merge"
            help="Merge the geometries in multiple JSON files into one."
            action=:command
    end

    parser["create"].description = "Create a geometry JSON file with obstructions of a specific type. Select a type for more help."

    for (as_string, constructor) in [
        ("walls", walls),
        ("cylinders", cylinders),
        ("spheres", spheres),
        ("annuli", annuli),
    ]
        add_arg_table!(parser["create"], as_string, Dict(
            :help => "Fill the geometry with $as_string.",
            :action => :command
        ))
        parser["create"][as_string].description = "Create a geometry JSON file filled with only $as_string with any properties defined by the flags."

        add_arg_table!(parser["create"][as_string],
            "number", Dict(
                :arg_type => Int,
                :help => "Number of obstructions to create.",
                :required => true,
            ),
            "output_file", Dict(
                :help => "Geometry JSON output filename.",
                :required => true,
            )
        )
        group = constructor(number=0)
        for unique_key in group.unique_keys
            field_value = group.field_values[unique_key]
            as_dict = Dict(
                :dest_name => String(unique_key),
                :help => field_to_docs(field_value),
                :required => field_value.field.required && isnothing(field_value.field.default_value),
                :default => field_value.field.default_value,
                :arg_type => field_type(field_value.field),
            )
            flag = "--" * String(unique_key)
            if field_type(field_value.field) == Bool
                for s in (:required, :default, :arg_type)
                    pop!(as_dict, s)
                end
                if field_value.field.default_value
                    flag = "--no-" * String(unique_key)
                    as_dict[:action] = :store_false
                else
                    as_dict[:action] = :store_true
                end
            elseif field_type(field_value.field) <: AbstractVector
                if field_value.field.only_group && field_type(field_value.field) <: SizedVector
                    as_dict[:nargs] = size(field_type(field_value.field))[1]
                else
                    as_dict[:nargs] = '+'
                    if !isnothing(as_dict[:default])
                        as_dict[:default] = [as_dict[:default]...]
                    end
                end
                as_dict[:arg_type] = eltype(field_type(field_value.field))
            elseif !field_value.field.only_group
                as_dict[:nargs] = '+'
                if !isnothing(as_dict[:default])
                    as_dict[:default] = [as_dict[:default]]
                end
            end
            add_arg_table!(parser["create"][as_string], flag, as_dict)
        end
    end

    add_arg_table!(parser["merge"],
        "output_file", Dict(
            :help => "A new geometry JSON file containing all the obstructions from the input files.",
            :required => true,
        ),
        "input_file", Dict(
            :help => "List of input geometry JSON files that should be merged. Each input file can contain obstructions of different types.",
            :required => true,
            :nargs => '+',
        ),
    )
    return parser
end

"""
    run_main([arguments])

Runs the `mcmr geometry` command line interface.
The supplied arguments can be provided as a sequence of strings (as provided by the terminal)
or as a dictionary (as provided by `ArgParse` after parsing).
By default it is set to `ARGS`.
"""
function run_main(args=ARGS::AbstractVector[<:AbstractString])
    parser = get_parser()
    return run_main(parse_args(args, parser))
end

function run_main(::Nothing)
    return Cint(1)
end

function run_main(args::Dict{<:AbstractString, <:Any})
    cmd = args["%COMMAND%"]
    if cmd == "create"
        run_create(args[cmd])
    end
    return Cint(0)
end

function run_create(args::Dict{<:AbstractString, <:Any})
    obstruction_type = args["%COMMAND%"]
    flags = args[obstruction_type]
    output_file = pop!(flags, "output_file")
    symbol_flags = Dict(Symbol(k) => v for (k, v) in flags)
    constructor = Dict(
        "walls" => walls,
        "cylinders" => cylinders,
        "spheres" => spheres,
        "annuli" => annuli,
    )[obstruction_type]
    result = constructor(;symbol_flags...)
    write_geometry(output_file, result)
end

end