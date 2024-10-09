"""
Types:
- [`GlobalProperties`](@ref)

Methods:
- [`R1`](@ref)
- [`R2`](@ref)
- [`off_resonance`](@ref)
- [`correct_for_timestep`](@ref)
- [`stick_probability`](@ref)
"""
module Properties

"""
    GlobalProperties(; R1=0, R2=0, off_resonance=0, permeability=0, dwell_time=0, surface_density=0, surface_relaxivity=0)

Stores global MRI and collision properties.
"""
struct GlobalProperties
    R1 :: Float64
    R2 :: Float64
    off_resonance :: Float64
end

function GlobalProperties(; R1=0, R2=0, off_resonance=0) 
    GlobalProperties(R1, R2, off_resonance)
end

for symbol in propertynames(GlobalProperties())
    @eval begin
        $symbol(global_properties::GlobalProperties) = getproperty(global_properties, $(QuoteNode(symbol)))
    end
end

function Base.show(io::IO, prop::GlobalProperties)
    print(io, "GlobalProperties(")
    for (func, units) in [
        (R1, "kHz"),
        (R2, "kHz"),
        (off_resonance, "kHz"),
    ]
        value = func(prop)
        if !iszero(value) && !isinf(value) && !isnan(value)
            print(io, "$(nameof(func))=$(value)$(units), ")
        end
    end
    print(io, ")")
end

"""
    correct_for_timestep(surface_relaxivity/permeability, timestep)

Corrects the surface relaxivity or permeability for the variability in the timestep during the simulation.

In Monte Carlo simulations the rate of collisions depends on the size of the timestep.
This means that as the timestep changes, the effect of surface relaxivity and permeability will depend on the timestep.
This function corrects the surface relaxivity and permeability values, so that their effect does not depend on timestep.
The user-provided surface relaxivity and permeability will be used unaltered if the timestep is 1 milliseconds.
"""
correct_for_timestep(probability, timestep) = 1 - (1 - probability)^sqrt(timestep)


"""
    stick_probability(surface_density, dwell_time, diffusivity, timestep)

Computes the probability of a spin getting stuck at the surface given
a [`surface_density`](@ref) and [`dwell_time`](@ref)
as well as the diffusivity (in um^2/ms) and the timestep (in ms).
"""
function stick_probability(surface_density::Number, dwell_time::Number, diffusivity::Number, timestep::Number)
    if iszero(surface_density)
        return 0.
    elseif iszero(dwell_time)
        error("Cannot have a dwell time of zero for any surface containing bound spins.")
    else
        return sqrt(π * timestep / diffusivity) * surface_density/ dwell_time / 2
    end
end
end