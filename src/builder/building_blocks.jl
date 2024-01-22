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
        if value isa GenericVariableRef || name == :builder
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


"""
    match_blocks!(block1, block2[, property_list])

Matches the listed properties between two [`BuildingBlock`](@ref) objects.
By default all shared properties (i.e., those with the same name) are matched.
"""
function match_blocks!(block1::BuildingBlock, block2::BuildingBlock, property_list)
    model = owner_model(block1)
    @assert model == owner_model(block2)
    for fn in property_list
        @constraint model fn(block1) == fn(block2)
    end
end

function match_blocks!(block1::BuildingBlock, block2::BuildingBlock)
    property_list = intersect(helper_functions(block1), helper_functions(block2))
    match_blocks!(block1, block2, property_list)
end


"""
Stores the parameters passed on a [`BuildingBlock`](@ref) constructor (of type `T`).

The parameters are temporarily stored in this format, until they can be added to a `SequenceBuilder`.

For example, the following
```julia
pg = PulsedGradient(qval=2)
```
will return a `BuildingBlockPlaceholder` object rather than a `PulsedGradient` object.
Only when this object is added to a `SequenceBuilder`, is the `PulsedGradient` actually initialised using the JuMP model of that `SequenceBuilder`:
```julia
sb = SequenceBuilder([pg])
```
You can access the initialised `PulsedGradient` through the `BuildingBlockPlaceholder` (`pg.concrete[]`) or directly through the `SequenceBuilder` (`sb[1]`)

Each Placeholder can be added to only a single `SequenceBuilder`, but it can be added multiple times to the same `SequenceBuilder`.
If added multiple times to the same `SequenceBuilder`, all variables will be matched between them.
"""
struct BuildingBlockPlaceholder{T<:BuildingBlock}
    args
    kwargs
    concrete :: Ref{T}
    BuildingBlockPlaceholder{T}(args...; kwargs...) where {T<:BuildingBlock} = new{T}(args, kwargs, Ref{T}())
end

function Base.show(io::IO, placeholder::BuildingBlockPlaceholder{T}) where {T}
    if isassigned(placeholder.concrete)
        print(io, "Assigned BuildingBlockPlaceholder for $(placeholder.concrete[])")
    else
        args = join(placeholder.args, ", ")
        kwargs = join(["$key=$value" for (key, value) in pairs(placeholder.kwargs)], ", ")
        print(io, "Unassigned BuildingBlockPlaceholder{$T}($args; $kwargs)")
    end
end


end