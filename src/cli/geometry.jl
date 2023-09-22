"""
Defines command line interface for `mcmr geometry`
"""
module Geometry

import ArgParse: ArgParseSettings, @add_arg_table!, add_arg_table!, parse_args, ArgParseError, usage_string
import StaticArrays: SizedVector
import ...Geometries.User.Obstructions: walls, Walls, spheres, Spheres, cylinders, Cylinders, annuli, Annuli, fields, field_to_docs, Field, ObstructionGroup, FieldValue
import ...Geometries.User.JSON: write_geometry, read_geometry

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
            elseif field_type(field_value.field) <: AbstractArray
                if field_value.field.only_group && field_type(field_value.field) <: SizedVector
                    as_dict[:nargs] = size(field_type(field_value.field))[1]
                else
                    as_dict[:nargs] = '+'
                    if !isnothing(as_dict[:default])
                        as_dict[:default] = [as_dict[:default]...]
                    end
                end
                as_dict[:arg_type] = eltype(field_type(field_value.field))
                if field_value.field.name == :rotation
                    pop!(as_dict, :arg_type)
                    as_dict[:default] = ["I"]
                end
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
    parse_user_argument(field_value, value, n_objects)

Parse the user argument to something that can actually be used internally by Julia.
"""
parse_user_argument(field_value::FieldValue{T}, value, n_objects) where {T} = value
function parse_user_argument(field_value::FieldValue{T}, value::Vector, n_objects) where {T}
    if field_value.field.name == :rotation
        return parse_user_rotation(value)
    end
    if length(value) == 0
        return nothing
    elseif T <: SizedVector
        nt = size(T)[1]
        if length(value) == nt
            return value
        elseif length(value) == nt * n_objects
            return [value[1 + i * nt: nt * (i + 1) for i in 0:n_objects-1]]
        else
            error("Expected $nt or $(nt * n_objects) values for $(field_value.name). Got $(length(value)) instead.")
        end
    else
        if length(value) == 1
            return value[1]
        elseif length(field_value) == n_objects
            return value
        end
    end
end


function parse_user_rotation(value::Vector)
    if length(value) == 1
        return Symbol(value[1])
    end
    parsed = parse.(Float64, value)
    if length(parsed) == 3
        return parsed
    else
        return reshape(parsed, (length(parsed) // 3, 3))
    end
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
    elseif cmd == "merge"
        run_merge(args[cmd])
    end
    return Cint(0)
end

function run_create(args::Dict{<:AbstractString, <:Any})
    obstruction_type = args["%COMMAND%"]
    flags = args[obstruction_type]
    output_file = pop!(flags, "output_file")
    constructor = Dict(
        "walls" => walls,
        "cylinders" => cylinders,
        "spheres" => spheres,
        "annuli" => annuli,
    )[obstruction_type]
    test_group = constructor(number=0)
    number = pop!(flags, "number")
    symbol_flags = Dict(Symbol(k) => parse_user_argument(getproperty(test_group, Symbol(k)), v, number) for (k, v) in flags)
    filtered = Dict(k=>v for (k, v) in symbol_flags if ~isnothing(v))
    result = constructor(;number=number, filtered...)
    write_geometry(output_file, result)
end

function run_merge(args)
    res = ObstructionGroup[]
    for input_file in args["input_file"]
        input_geom = read_geometry(input_file)
        if input_geom isa ObstructionGroup
            push!(res, input_geom)
        else
            append!(res, input_geom)
        end
    end
    write_geometry(args["output_file"], res)
end

end