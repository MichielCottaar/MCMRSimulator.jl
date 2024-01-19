module BuildingBlocks
import JuMP: has_values, GenericVariableRef, value, Model, @constraint, @objective, owner_model
import ...Sequences: RFPulse, InstantRFPulse, MRGradients, InstantGradient, Sequence
import ...Scanners: Scanner

"""
Parent type for all individual components out of which a sequence can be built.

Required methods:
- [`duration`](@ref)(block, parameters): returns block duration in ms
- [`scanner_constraints!`](@ref)(model, block, scanner): adds scanner constraints to the block
- [`to_components`](@ref)(block, parameters): converts the block into components recognised by the MCMR simulator
- [`helper_functions`](@ref): A list of all functions that are used to compute properties of the building block. Any of these can be used in constraints or objective functions.
"""
abstract type BuildingBlock end

function owner_model(bb::BuildingBlock)
    owner_model(duration(bb))
end

has_values(block::BuildingBlock) = has_values(owner_model(block))


"""
    duration(building_block)

The duration of the building block in ms.
"""
function duration end

"""
    to_components(building_block, parameters)

Converts the building block into components recognised by the MCMR simulator. These components are:
- [`RFPulse`](@ref)
- [`InstantRFPulse`](@ref)
- [`MRGradients`](@ref)
- [`InstantGradient`](@ref)
"""
function to_components end

"""
    scanner_constraints!([model, ]building_block, scanner)

Adds any constraints from a specific scanner to a sequence or BuildingBlock
"""
function scanner_constraints!(building_block::BuildingBlock, scanner::Scanner)
    scanner_constraints!(owner_model(building_block), building_block, scanner)
end

"""
    helper_functions(building_block)

Returns a list of function that can be called to constrain the `building_block`.
"""
helper_functions(bb::BuildingBlock) = helper_functions(typeof(bb))


struct _BuildingBlockPrinter
    bb :: BuildingBlock
    number :: Integer
end

function Base.show(io::IO, block::BuildingBlock)
    print(io, string(typeof(block)), "(")
    for name in propertynames(block)
        value = getproperty(block, name)
        if value isa GenericVariableRef
            continue
        end
        print(io, name, "=", repr(value), ", ")
    end

    if has_values(block)
        for fn in helper_functions(block)
            print(io, "$(nameof(fn))=$(value(fn(block))), ")
        end
    end
    print(io, ")")
end


"""
    set_simple_constraints!(model, block)

Add any constraints or objective functions to the properties of a [`BuildingBlock`](@ref) in the JuMP `model`.

Each keyword argument has to match one of the functions in [`helper_functions`](@ref)(block).
If set to a numeric value, a constraint will be added to fix the function value to that numeric value.
If set to `:min` or `:max`, minimising or maximising this function will be added to the cost function.
"""
function set_simple_constraints!(model::Model, block::BuildingBlock, kwargs)
    to_funcs = Dict(nameof(fn) => fn for fn in helper_functions(block))
    for (key, value) in kwargs
        apply_simple_constraint!(model, to_funcs[key](block), value)
    end
    nothing
end

"""
    apply_simple_constraint!(model, variable, value)

Add a single constraint or objective to the JuMP `model`.
This is an internal function used by [`set_simple_constraints`](@ref).
"""
apply_simple_constraint!(model::Model, variable, ::Val{:min}) = @objective model Min objective_function(model) + variable
apply_simple_constraint!(model::Model, variable, ::Val{:max}) = @objective model Min objective_function(model) - variable
apply_simple_constraint!(model::Model, variable, value::Number) = @constraint model variable == value


end