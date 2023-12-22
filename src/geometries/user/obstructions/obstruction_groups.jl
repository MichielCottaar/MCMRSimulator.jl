module ObstructionGroups

import ..ObstructionTypes: ObstructionType, fields, key_value_pairs
import ..Fields: Field, FieldValue, isglobal


struct ObstructionGroup{N}
    type :: ObstructionType{N}
    field_values :: Dict
    unique_keys :: Vector{Symbol}
    n_obstructions :: Int
end


Base.length(g::ObstructionGroup) = g.n_obstructions
function Base.getindex(g::ObstructionGroup, index::Int)
    return IndexedObstruction(g, index)
end
Base.iterate(g::ObstructionGroup) = iterate(g, 1)
Base.iterate(g::ObstructionGroup, state::Int) = state > g.n_obstructions ? nothing : (g[state], state + 1)

struct IndexedObstruction{N}
    group :: ObstructionGroup{N}
    index :: Int
end


Base.propertynames(og::ObstructionGroup) = og.unique_keys
Base.propertynames(io::IndexedObstruction) = propertynames(io.group)
function Base.getproperty(og::ObstructionGroup, s::Symbol) 
    d = getfield(og, :field_values)
    if s in keys(d)
        return d[s]
    else
        return getfield(og, s)
    end
end

function Base.setproperty!(og::ObstructionGroup, s::Symbol, value) 
    og.field_values[s].value = value
end

function Base.getproperty(io::IndexedObstruction, s::Symbol) 
    if s in (:group, :index)
        return getfield(io, s)
    end
    return io.group.field_values[s][io.index]
end

function Base.setproperty!(io::IndexedObstruction, s::Symbol, value) 
    io.group.field_values[s][io.index] = value
end

function ObstructionGroup(type::ObstructionType; number=nothing, kwargs...)

    (unique_keys, field_values) = key_value_pairs(type, isnothing(number) ? 0 : number)
    for (key, value) in kwargs
        field_values[key].value = value
    end

    if isnothing(number)
        # try to guess the number of obstructions
        for fv in values(field_values)
            if ~isglobal(fv)
                number = length(fv.value)
                println("There are $number $(type) based on values for $(fv.field)")
                break
            end
        end
        if isnothing(number)
            number = 1
        end
        for fv in values(field_values)
            fv.n_obstructions = number
        end
    end

    group = ObstructionGroup(type, field_values, unique_keys, number)
    return group
end

function Base.show(io::IO, og::ObstructionGroup)
    print(io, "$(og.n_obstructions) $(lowercase(String(og.type.plural)))(")
    not_set = Symbol[]
    required_missing = Symbol[]
    for key in propertynames(og)
        fv = og.field_values[key]
        if isnothing(fv.value)
            if fv.field.required
                push!(required_missing, key)
            else
                push!(not_set, key)
            end
        elseif (isglobal(fv) && key != :vertices) || og.n_obstructions < 5
            print(io, "$key = $(fv.value), ")
        else
            print(io, "$key = <unique values>, ")
        end
    end
    print(io, ")")
    if length(required_missing) > 0
        text = join(required_missing, ", ", " and ")
        print(io, "\nMissing required fields: $text")
    end
    if length(not_set) > 0
        text = join(not_set, ", ", " and ")
        print(io, "\nNot set: $text")
    end
end

function Base.show(io::IO, obstruction::IndexedObstruction)
    og = obstruction.group
    print(io, "$(og.type.singular) at index $(obstruction.index) (")
    already_done = Set()
    for (key, fv) in og.field_values
        if fv in already_done
            continue
        end
        push!(already_done, fv)
        my_value = fv[obstruction.index]
        if ~isnothing(my_value)
            print(io, "$key = $(my_value), ")
        end
    end
    print(io, ")")
end

end