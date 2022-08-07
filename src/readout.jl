mutable struct Simulation{N, M, T<:AbstractFloat}
    # N sequences, M spins, datatype T
    sequences :: SVector{N, <:Sequence}
    micro::Microstructure
    timestep::T
    regular :: AbstractVector{MultiSnapshot{N, M, T}}
    store_every :: T
    readout :: SVector{N, AbstractVector{Snapshot{M, T}}}
    latest :: MultiSnapshot
    function Simulation(spins, sequences :: AbstractVector{<:Sequence}, micro::Microstructure; store_every :: Real=5., timestep :: Real=0.5)
        if isa(spins, Spin)
            spins = [spins]
        end
        if isa(spins, AbstractVector{<:Spin})
            snap = Snapshot(spins, 0.)
        end
        if isa(spins, Snapshot)
            snap = MultiSnapshot(snap, nseq)
        end
        nseq = length(sequences)
        nspins = length(snap.spins)
        snap = MultiSnapshot(Snapshot(spins, 0.), nseq)
        new{nseq, nspins, typeof(timestep)}(
            SVector{nseq}(sequences),
            micro,
            timestep,
            MultiSnapshot{nseq, nspins, Float64}[],
            store_every,
            [Snapshot{nspins, Float64}[] for _ in 1:nseq],
            snap
        )
    end
end