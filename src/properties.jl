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

for (symbol, summary, opposite) in [
    (:R1, "longitudinal relaxation rate (kHz)", "T1"),
    (:T1, "longitudinal relaxation time (ms)", "R1"),
    (:R2, "transverse relaxation rate (kHz)", "T2"),
    (:T2, "transverse relaxation time (ms)", "R2"),
    (:off_resonance, "off-resonance field (kHz)", nothing)
]
    symbol_str = string(symbol)

    exclusive_statement = isnothing(opposite) ? "" : """ It is the inverse of $(opposite). $(opposite) and $(symbol_str) cannot both be set."""
    @eval begin
        """
            $($symbol_str)(properties)

        The $($summary) stored in [`MRIProperties`](@ref).$($exclusive_statement)

        It can be set by the user when creating a [`Simulation`](@ref), in which case it will be stored in the [`GlobalProperties`](@ref).
        When creating objects in this geometry, the $($symbol_str) can be updated in two ways:
        - setting the `$($symbol_str)_inside` keyword affecting spins within the object
        - setting the `$($symbol_str)_surface` keyword affecting spins stuck on the surface of the object
        Both will be stored within a [`ObstructionProperties`](@ref) object.
        """
        function $(symbol) end
    end
end

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
    CollisionProperties(MT_fraction, permeability, surface_density, dwell_time)

Object used to define the collision properties within a [`ObstructionProperties`](@ref) or a [`GlobalProperties`] object.
Values of NaN are used in [`ObstructionProperties`] to indicate that the default collision parameter should be used.

These properties can be retrieved using:
- [`MT_fraction`](@ref)
- [`permeability`](@ref)
- [`surface_density`](@ref)
- [`dwell_time`](@ref)
"""
struct CollisionProperties
    MT_fraction :: Float
    permeability :: Float
    surface_density :: Float
    dwell_time :: Float
end

"""
    ObstructionProperties(; MT_fraction, permeability, relative_density, dwell_time, R1_inside, T1_inside, R2_inside, T2_inside, off_resonance_inside, R1_surface, T1_surface, R2_surface, T2_surface, off_resonance_surface)

Defines the collision and relaxation properties for an obstruction.
Either the relaxation rate (R1/R2) or the relaxation time (T1/T2) should be set (or neither).

Collision parameters can be retreived using their accessors: [`MT_fraction`](@ref), [`permeability`](@ref), [`surface_density`](@ref), and [`dwell_time`](@ref).

MRI relaxation parameters within the obstruction can be retrieved by calling their accessor on `obstruction_properties.inside` (for spins inside the obstruction) or `obstruction_properties.surface` (for spins stuck at the surface of the obstruction).
The accessors are [`T1`](@ref), [`R1`](@ref), [`T2`](@ref), [`R2`](@ref), and [`off_resonance`](@ref).

For any property not set (or set to NaN) the default values will be used (i.e., from the [`GlobalProperties`](@ref), which is set during the creation of the [`Simulation`](@ref)).
"""
struct ObstructionProperties
    # In future add MRIProperties for particles stuck on the surface
    collision :: CollisionProperties
    inside :: MRIProperties
    surface :: MRIProperties
    id :: UUID
end

function ObstructionProperties(;
    MT_fraction::Number=NaN,
    permeability::Number=NaN,
    R1_inside::Number=NaN, T1_inside::Number=NaN,
    R2_inside::Number=NaN, T2_inside::Number=NaN,
    off_resonance_inside::Number=NaN,
    R1_surface::Number=NaN, T1_surface::Number=NaN,
    R2_surface::Number=NaN, T2_surface::Number=NaN,
    off_resonance_surface::Number=NaN,
    surface_density::Number=NaN,
    dwell_time::Number=NaN,
)
    ObstructionProperties(
        CollisionProperties(Float(MT_fraction), Float(permeability), Float(surface_density), Float(dwell_time)),
        MRIProperties(; R1=R1_inside, T1=T1_inside, R2=R2_inside, T2=T2_inside, off_resonance=off_resonance_inside),
        MRIProperties(; R1=R1_surface, T1=T1_surface, R2=R2_surface, T2=T2_surface, off_resonance=off_resonance_surface),
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
- [`surface_density`](@ref)
- [`dwell_time`](@ref)
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
        surface_density::Number=0,
        dwell_time::Number=NaN,
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
            CollisionProperties(Float(MT_fraction), Float(permeability), Float(surface_density), Float(dwell_time)),
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
        (surface_density, "um"),
        (dwell_time, "ms"),
    ]
        value = func(prop)
        if !iszero(value) && !isinf(value) && !isnan(value)
            print(io, "$(nameof(func))=$(value)$(units), ")
        end
    end
    print(io, ")")
end

for symbol in (:permeability, :MT_fraction, :surface_density, :dwell_time)
    @eval $(symbol)(c::CollisionProperties) = c.$(symbol)
    @eval $(symbol)(o::ObstructionProperties) = $(symbol)(o.collision)
    @eval $(symbol)(g::GlobalProperties) = $(symbol)(g.collision)

    @eval begin
        """
            $($symbol)(properties)
            $($symbol)(obstruction_properties, global_properties)

        The $($symbol), which affects the spin when colliding with an obstruction.
        It can be either set when creating a [`Simulation`](@ref), in which case it is stored in [`GlobalProperties`](@ref)
        or when creating any [`Obstruction`](@ref), in which case it is stored in [`ObstructionProperties`](@ref).
        """
        function $(symbol)(o::ObstructionProperties, defaults) 
            value = $(symbol)(o)
            if isnan(value)
                value = $(symbol)(defaults)
            end
            @assert !isnan(value)
            return value
        end
    end
end

for symbol in (:R1, :T1, :R2, :T2, :off_resonance)
    @eval $(symbol)(g::GlobalProperties) = $(symbol)(g.mri)
end

"""
    stick_probability(surface_density, dwell_time, diffusivity, timestep)
    stick_probability(properties, diffusivity, timestep)

Computes the probability of a spin getting stuck at the surface given
a [`surface_density`](@ref) and [`dwell_time`](@ref) from a [`CollisionProperties`](@ref)
as well as the diffusivity (in um^2/ms) and the timestep (in ms).
"""
function stick_probability(surface_density::Number, dwell_time::Number, diffusivity::Number, timestep::Number)
    return sqrt(π * timestep / diffusivity) * surface_density/ dwell_time / 2
end

function stick_probability(properties, diffusivity::Number, timestep::Number)
    return stick_probability(surface_density(properties), dwell_time(properties), diffusivity, timestep)
end

"""
    correct_for_timestep(MT_fraction/permeability, timestep)

Corrects the MT fraction or permeability for the variability in the timestep during the simulation.

In Monte Carlo simulations the rate of collisions depends on the size of the timestep.
This means that as the timestep changes, the effect of MT fraction and permeability will depend on the timestep.
This function corrects the MT fraction and permeability values, so that their effect does not depend on timestep.
The user-provided MT fraction and permeability will be used as is if the timestep is 1 milliseconds.
"""
correct_for_timestep(probability, timestep) = 1 - (1 - probability)^sqrt(timestep)