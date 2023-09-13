"""
Defines command line interface for `mcmr geometry`
"""
module Geometry

import ArgParse: ArgParseSettings, @add_arg_table!, add_arg_table!, parse_args
import ...Geometries.User.Obstructions: walls, Walls, spheres, Spheres, cylinders, Cylinders, annuli, Annuli, fields, field_to_docs, Field

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
            required = field_value.field.required & isnothing(field_value.field.default_value)
            if field_type(field_value.field) == Bool
                if field_value.field.default_value
                    add_arg_table!(parser["create"][as_string], "--no-" * String(unique_key), Dict(
                        :action => :store_false,
                        :dest_name => String(unique_key),
                        :help => field_to_docs(field_value)
                    ))
                else
                    add_arg_table!(parser["create"][as_string], "--" * String(unique_key), Dict(
                        :action => :store_true,
                        :help => field_to_docs(field_value)
                    ))
                end
            else
                add_arg_table!(parser["create"][as_string], "--" * String(unique_key), Dict(
                    :help => field_to_docs(field_value),
                    :required => required,
                    :default => field_value.field.default_value
                ))
            end
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
    run_main(parse_args(args, parser))
end

function run_main(::Nothing)
    return
end

function run_main(args::Dict{<:AbstractString, <:Any})
    cmd = args["%COMMAND%"]
    @show args[cmd]
end


end