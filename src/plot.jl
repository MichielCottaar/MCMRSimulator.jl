
@userplot SpinQuiver
@recipe function f(sq::SpinQuiver)
    res = sq.args[1]
    @assert isa(res, AbstractVector{Spin})
    coords = [map(x -> x.position[1:2], res)]
    seriestype := :quiver
    #z := map(x -> x.position[3], res)
    quiver := Tuple([map(x -> x[index], vector.(res)) for index in 1:2])
    aspect_ratio --> :equal
    map(x -> x.position[1], res), map(x -> x.position[2], res)
end
