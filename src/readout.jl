mutable struct Simulation{N, T<:AbstractFloat, M<:Microstructure}
    # N sequences, datatype T
    sequences :: SVector{N, <:Sequence}
    micro::M
    timestep::T
    regular :: AbstractVector{MultiSnapshot{N, T}}
    store_every :: T
    readout :: SVector{N, AbstractVector{Snapshot{T}}}
    latest :: MultiSnapshot{N, T}
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
        snap = MultiSnapshot(Snapshot(spins, 0.), nseq)
        new{nseq, typeof(timestep), typeof(micro)}(
            SVector{nseq}(sequences),
            micro,
            timestep,
            MultiSnapshot{nseq, Float64}[],
            store_every,
            [Snapshot{Float64}[] for _ in 1:nseq],
            snap
        )
    end
end