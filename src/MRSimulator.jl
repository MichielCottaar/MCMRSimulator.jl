module MRSimulator
import StaticArrays: SA_F64, MVector, SVector
using LinearAlgebra
import Base
import RecipesBase: RecipesBase, @userplot, @recipe


include("constants.jl")
include("spin.jl")
include("field.jl")
include("sequence.jl")
include("plot.jl")

function evolve_to_time(spin :: Spin, micro :: Microstructure, new_time :: Real, B0=3.)
    if spin.time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time > spin.time
        timestep = new_time - spin.time
        spin = Spin(
            new_time,
            spin.position,
            relax(spin.orientation, micro(spin.position), timestep, B0)
        )
    end
    spin
end

function evolve(spin :: Spin, micro :: Microstructure, sequence :: Sequence; store_every=1., B0=3., nTR=1)
    sequence_index = 1
    readout_index = 1
    spins = typeof(spin)[]
    times = MVector{3, Float64}([
        nextfloat(nTR * sequence.TR),
        time(sequence, sequence_index),
        (readout_index - 1) * store_every
    ])
    while true
        next = argmin(times)
        spin = evolve_to_time(spin, micro, times[next], B0)
        if next == 1
            break
        elseif next == 2
            spin = apply(sequence[sequence_index], spin)
            sequence_index += 1
            times[2] = time(sequence, sequence_index)
        elseif next == 3
            push!(spins, spin)
            readout_index += 1
            times[3] = (readout_index - 1) * store_every
        end
    end
    spins
end

function evolve(spins :: Vector{Spin}, micro :: Microstructure, sequence :: Sequence; store_every=1., B0=3., nTR=1)
    nspins = length(spins)
    result = fill(Spin[], nspins)
    Threads.@threads for i = 1:nspins
        result[i] = evolve(spins[i], micro, sequence, store_every=store_every, B0=B0, nTR=nTR)
    end
    hcat(result...)
end

export evolve, Microstructure, field, Sequence, RFPulse, Spin, vector
end

