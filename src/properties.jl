"""
    MRIProperties(; R1=, R2=, T1=, T2=, off_resonance=)

Object used to define the MRI relaxation properties within a [`ObstructionProperties`](@ref) or a [`GlobalProperties`] object.
Either the relaxation rate (R1/R2) or the relaxation time (T1/T2) should be set (or neither).

Values of NaN are used in [`ObstructionProperties`] to indicate that the default relaxation parameter should be used.

These relaxation variables can be retrieved using:
- [`T1`](@ref)
- [`R1`](@ref)
- [`T2`](@ref)
- [`R2`](@ref)
- [`off_resonance`](@ref)
"""
struct MRIProperties
    R1 :: Float
    R2 :: Float
    off_resonance :: Float
end

function MRIProperties(; R1=NaN, R2=NaN, T1=NaN, T2=NaN, off_resonance=NaN)
    if isnan(R1)
        R1 = isnan(T1) ? NaN : 1/T1
    elseif !isnan(T1)
        error("R1 and T1 cannot both be set at the same time")
    end
    if isnan(R2)
        R2 = isnan(T2) ? NaN : 1/T2
    elseif !isnan(T2)
        error("R2 and T2 cannot both be set at the same time")
    end
    MRIProperties(Float(R1), Float(R2), Float(off_resonance))
end

R1(p::MRIProperties) = p.R1
R2(p::MRIProperties) = p.R2
T1(p::MRIProperties) = 1/R1(p)
T2(p::MRIProperties) = 1/R2(p)
off_resonance(p::MRIProperties) = p.off_resonance

"""
    empty_mri_properties(properties)

Returns true if none of the parameters have been set
"""
empty_mri_properties(p::MRIProperties) = isnan(R1(p)) && isnan(R2(p)) && isnan(off_resonance(p))


"""
    merge_mri_parameters(properties[, default_values=<empty>])

Merges multiple MRI properties into one.
An error is raised if they are inconsistent with each other
"""
function merge_mri_parameters(properties, default_values=nothing)
    if isnothing(default_values)
        R1_val = R2_val = off_resonance_val = NaN
    else
        R1_val = default_values.R1
        R2_val = default_values.R2
        off_resonance_val = default_values.off_resonance
    end
    set_R1 = false
    set_R2 = false
    set_off_resonance = false
    for prop in properties
        if !isnan(prop.R1)
            if set_R1 && !(R1_val ≈ R1(prop))
                error("Inconsistent values found for R1 when mergine MRI properties")
            end
            R1_val = R1(prop)
            set_R1 = true
        end
        if !isnan(prop.R2)
            if set_R2 && !(R2_val ≈ R2(prop))
                error("Inconsistent values found for R2 when mergine MRI properties")
            end
            R2_val = R2(prop)
            set_R2 = true
        end
        if !isnan(prop.off_resonance)
            if set_off_resonance && !(off_resonance_val ≈ off_resonance(prop))
                error("Inconsistent values found for off_resonance when mergine MRI properties")
            end
            off_resonance_val = off_resonance(prop)
            set_off_resonance = true
        end
    end
    return MRIProperties(R1_val, R2_val, off_resonance_val)
end

function setproperty!(props::MRIProperties, symbol, value)
    if symbol == :T1
        setfield!(props, :R1, 1/value)
    end
    if symbol == :T2
        setfield!(props, :R2, 1/value)
    end
    setfield!(props, symbol, value)
end



"""
    CollisionProperties(MT_fraction, permeability)

Object used to define the collision properties within a [`ObstructionProperties`](@ref) or a [`GlobalProperties`] object.
Values of NaN are used in [`ObstructionProperties`] to indicate that the default collision parameter should be used.

These properties can be retrieved using:
- [`MT_fraction`](@ref)
- [`permeability`](@ref)
"""
struct CollisionProperties
    MT_fraction :: Float
    permeability :: Float
end

"""
    ObstructionProperties(; MT_fraction=NaN, permeability=NaN, R1_inside=NaN, T1_inside=NaN, R2_inside=NaN, T2_inside=NaN, off_resonance_inside=NaN)

Defines the collision and relaxation properties for an obstruction.
Either the relaxation rate (R1/R2) or the relaxation time (T1/T2) should be set (or neither).

Collision parameters can be retreived using their accessors: [`MT_fraction`](@ref), and [`permeability`](@ref).

MRI relaxation parameters within the obstruction can be retrieved by calling their accessor on `obstruction_properties.inside`.
The accessors are [`T1`](@ref), [`R1`](@ref), [`T2`](@ref), [`R2`](@ref), and [`off_resonance`](@ref).

Values of NaN are used internally to indicate that the default parameters (from [`GlobalProperties`](@ref)) should be used.
"""
struct ObstructionProperties
    # In future add MRIProperties for particles stuck on the surface
    collision :: CollisionProperties
    inside :: MRIProperties
    id :: UUID
end

function ObstructionProperties(;
    MT_fraction::Number=NaN,
    permeability::Number=NaN,
    R1_inside::Number=NaN, T1_inside::Number=NaN,
    R2_inside::Number=NaN, T2_inside::Number=NaN,
    off_resonance_inside::Number=NaN,
)
    ObstructionProperties(
        CollisionProperties(Float(MT_fraction), Float(permeability)),
        MRIProperties(; R1=R1_inside, T1=T1_inside, R2=R2_inside, T2=T2_inside, off_resonance=off_resonance_inside),
        uuid1()
    )
end

"""
    GlobalProperties(; MT_fraction=0, permeability=0, R1=0, T1=Inf, R2=0, T2=Inf, off_resonance=0)

Defines a default set of collision and relaxation parameters for a simulation.
These parameters can be locally overriden by [`ObstructionProperties`](@ref).

Parameters can be accessed using their accessors:
- [`MT_fraction`](@ref)
- [`permeability`](@ref)
- [`T1`](@ref)
- [`R1`](@ref)
- [`T2`](@ref)
- [`R2`](@ref)
- [`off_resonance`](@ref)
"""
struct GlobalProperties
    collision :: CollisionProperties
    mri :: MRIProperties
    function GlobalProperties(;
        MT_fraction::Number=0,
        permeability::Number=0,
        R1::Number=NaN, T1::Number=NaN,
        R2::Number=NaN, T2::Number=NaN,
        off_resonance::Number=0,
    )
        if isnan(R1) & isnan(T1)
            R1 = 0
        end
        if isnan(R2) & isnan(T2)
            R2 = 0
        end
        res = new(
            CollisionProperties(Float(MT_fraction), Float(permeability)),
            MRIProperties(; R1=R1, T1=T1, R2=R2, T2=T2, off_resonance=off_resonance),
        )
        return res
    end
end

function Base.show(io::IO, prop::GlobalProperties)
    print(io, "GlobalProperties(")
    for (func, units) in [
        (T1, "ms"),
        (T2, "ms"),
        (off_resonance, "kHz"),
        (MT_fraction, ""),
        (permeability, ""),
    ]
        value = func(prop)
        if !iszero(value) && !isinf(value)
            print(io, "$(nameof(func))=$(value)$(units), ")
        end
    end
    print(io, ")")
end

for symbol in (:permeability, :MT_fraction)
    @eval $(symbol)(c::CollisionProperties) = c.$(symbol)
    @eval $(symbol)(o::ObstructionProperties) = $(symbol)(o.collision)
    @eval $(symbol)(g::GlobalProperties) = $(symbol)(g.collision)
end

for symbol in (:R1, :T1, :R2, :T2, :off_resonance)
    @eval $(symbol)(g::GlobalProperties) = $(symbol)(g.mri)
end