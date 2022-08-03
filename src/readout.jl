struct TR_Readout
    sequence :: Sequence
    regular :: AbstractVector{Snapshot}
    store_every :: Real
    data :: AbstractVector{Snapshot}
    final :: Snapshot
    all :: AbstractVector{Snapshot}
    function TR_Readout(sequence :: Sequence, nspins :: Int, store_every :: Real)
        nregular = iszero(store_every) ? 0 : Int(div(sequence.TR, store_every, RoundDown))
        regular = [Snapshot(nspins, (index - 1) * store_every) for index in 1:nregular]
        data = [Snapshot(nspins, p.time) for p in sequence.pulses if isa(p, Readout)]
        final = Snapshot(nspins, sequence.TR)
        new(
            sequence,
            regular,
            store_every,
            data,
            final, 
            sort!([data..., regular..., final], by=s->time(s))
        )
    end
end

struct SpinReadout
    parent :: TR_Readout
    index :: Integer
end


Base.getindex(ro::TR_Readout, i) = ro.all[i]
Base.iterate(ro::TR_Readout, state :: Int) = state > length(ro) ? nothing : (ro.all[state], state + 1)
Base.iterate(ro::TR_Readout) = iterate(ro, 1)
Base.eltype(::Type{TR_Readout}) = Snapshot
for param in (:length, :lastindex)
    @eval Base.$param(ro::TR_Readout) = $param(ro.all)
end

Base.getindex(ro::SpinReadout, i::Int) = TR_Readout.all[i][ro.index]
function Base.setindex!(ro::SpinReadout, value::Spin, i::Int) 
    ro.parent.all[i].spins[ro.index] = value
end
