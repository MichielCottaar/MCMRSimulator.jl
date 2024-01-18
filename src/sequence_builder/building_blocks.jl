module BuildingBlocks
import JuMP: has_values, GenericVariableRef, value, Model, @constraint, @objective
import ...Sequences: RFPulse, InstantRFPulse, MRGradients, InstantGradient, Sequence
import ...Scanners: Scanner

"""
Parent type for all individual components out of which a sequence can be built.

Required methods:
- [`duration`](@ref): returns block duration in ms
- [`to_components`](@ref): converts the block into components recognised by the MCMR simulator
- [`apply_scanner!`](@ref): apply [`Scanner`](@ref) constraints to the building block.
- [`helper_functions`](@ref): A list of all functions that are used to compute properties of the building block. Any of these can be used in constraints or objective functions.
- [`show_helper_functions`](@ref): A list of all functions and unit names that should be used when printing the building block.
"""
abstract type BuildingBlock end

"""
    duration(building_block)

The duration of the building block in ms.
"""
function duration end

"""
    to_components(building_block)

Converts the building block into components recognised by the MCMR simulator. These components are:
- [`RFPulse`](@ref)
- [`InstantRFPulse`](@ref)
- [`MRGradients`](@ref)
- [`InstantGradient`](@ref)
"""
function to_components end


"""
    apply_scanner!([model, ]scanner, block)

Apply any constraints from a specific MRI [`Scanner`](@ref) to the properties of the [`BuildingBlock`](@ref).
"""
function apply_scanner!(scanner::Scanner, block::BuildingBlock)
    apply_scanner!(get_model(block), scanner, block)
end


"""
    helper_functions(building_block)

Returns a list of function that can be called to constrain the `building_block`.
"""
helper_functions(bb::BuildingBlock) = helper_functions(typeof(bb))

"""
    show_helper_functions(building_block)

Returns a list of functions and user_names that are used to print the state of the `building_block`.
"""
show_helper_functions(bb::BuildingBlock) = show_helper_functions(typeof(bb))



"""
Raised by [`get_model`](@ref) if the JuMP model can not be found for a specific building block.
"""
struct NoModelException <: Exception end

function get_model(block::BuildingBlock)
    for param_name in propertynames(block)
        param = getproperty(block, param_name)
        if param isa GenericVariableRef
            return param.model
        elseif param isa BuildingBlock
            try
                return get_model(param)
            catch e
                if isa(e, NoModelException)
                    continue
                else
                    throw(e)
                end
            end
        end
    end
    throw(NoModelException, "Could not find underlying JuMP model for $(typeof(comp))")
end

function has_values(block::BuildingBlock)
    try
        return has_values(get_model(block))
    catch e
        if e isa NoModelException
            return true
        else
            throw(e)
        end
    end
end

function Base.show(io::IO, block::BuildingBlock)
    print(io, typeof(block))
    if has_values(block)
        print(io, "(")
        for (func, units) in show_helper_functions(block)
            print(io, nameof(func))
            print(io, "=")
            print(io, value(func(block)))
            print(io, units)
            print(io, ", ")
        end
        print(io, ")")
    end
end


"""
    set_kwargs_constraints!(model, block; kwargs...)

Add any constraints or objective functions to the properties of a [`BuildingBlock`](@ref) in the JuMP `model`.

Each keyword argument has to match one of the functions in [`helper_functions`](@ref)(block).
If set to a numeric value, a constraint will be added to fix the function value to that numeric value.
If set to `:min` or `:max`, minimising or maximising this function will be added to the cost function.
"""
function set_kwargs_constraints!(model::Model, block::BuildingBlock; kwargs...)
    to_funcs = Dict(nameof(func) => func for func in helper_functions(typeof(block)))
    for key in keys(kwargs)
        if !(key in keys(to_funcs))
            error("unrecognised parameter $key passed on PulsedGradient constructor")
        end
        value = kwargs[key]
        func = to_funcs[key]
        apply_constraint!(model, block, to_funcs[key], kwargs[key])
    end
    nothing
end

"""
    apply_constraint!(model, building_block, func, value)

Add a single constraint or objective to the JuMP `model`.
This is an internal function used by [`set_kwargs_constraints`](@ref).
"""
apply_constraint!(model::Model, block::BuildingBlock, func::Function, value::Nothing) = nothing
apply_constraint!(model::Model, block::BuildingBlock, func::Function, value::Symbol) = apply_constraint!(model, func, block, Val(value))
apply_constraint!(model::Model, block::BuildingBlock, func::Function, value::Val{:min}) = @objective model Min objective_function(model) + func(block)
apply_constraint!(model::Model, block::BuildingBlock, func::Function, value::Val{:max}) = @objective model Min objective_function(model) - func(block)
apply_constraint!(model::Model, block::BuildingBlock, func::Function, value::Number) = @constraint model func(block) == value


"""
    add_parameters_to_docs(building_block)

Return updated documentation of a [`BuildingBlock`](@ref) type that includes descriptions of all the helper functions in [`helper_functions`](@ref).
"""
function add_parameters_to_docs(building_block::Type{<:BuildingBlock}, initial_docs)
    for func in helper_functions(building_block)
        new_docs = string(Base.doc(func))
        stripped_docs = join(split(new_docs, '\n')[5:end], ' ')
        initial_docs = initial_docs * "- [`$(nameof(func))`](@ref): $(stripped_docs)\n"
    end
    return initial_docs
end


end