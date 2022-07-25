
@recipe function _f(snap::Snapshot)
    res = snap.spins
    @assert isa(res, AbstractVector{Spin})
    coords = [map(x -> x.position[1:2], res)]
    seriestype := :quiver
    #z := map(x -> x.position[3], res)
    quiver := Tuple([map(x -> x[index], vector.(res)) for index in 1:2])
    aspect_ratio --> :equal
    map(x -> x.position[1], res), map(x -> x.position[2], res)
end


@recipe function _f(sequence::Sequence)
    @assert isa(sequence, Sequence)
    total_time = sequence.TR
    color --> "black"
    times = map(t -> t.time, sequence.pulses)
    height = map(t -> flip_angle(t), sequence.pulses)
    max_height = maximum(height)
    label := nothing
    @series begin
        seriestype := :sticks
        lw --> 5
        times, height
    end
    @series begin
        seriestype := :scatter
        shape := :utriangle
        msize --> 10
        times, height
    end
    @series begin
        seriestype := :scatter
        alpha := 0
        times, height .+ max_height / 10
    end
    seriestype := :scatter
    annotations := [(p.time, flip_angle(p) + max_height / 10, string(Int(round(flip_angle(p))))) for p in sequence.pulses]
    ([], [])
end