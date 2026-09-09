"""Allocation-free views over cached hierarchical inside indices."""
module InsideViews

export InsideView, child_view, has_other

struct InsideView{V<:AbstractVector}
    indices::V
    first::Int
    last::Int
    depth::Int
end

InsideView(indices::V) where {V<:AbstractVector} = InsideView(indices, 1, length(indices), 0)

Base.length(view::InsideView) = max(view.last - view.first + 1, 0)
Base.isempty(view::InsideView) = view.first > view.last

function _local_component(view::InsideView, index, level)
    index[view.depth + level]
end

function Base.iterate(view::InsideView, state=view.first)
    state > view.last && return nothing
    index = view.indices[state]
    local_index = view.depth == 0 ? index : index[(view.depth + 1):end]
    local_index, state + 1
end

function _matches(view::InsideView, index, prefix)
    length(index) >= view.depth + length(prefix) || return false
    all(_local_component(view, index, level) == prefix[level] for level in eachindex(prefix))
end

child_view(::Val{false}, _, _) = nothing

child_view(::Val{true}, ::Nothing, ::Tuple) = nothing

function child_view(::Val{true}, view::InsideView, prefix::Tuple)
    isempty(prefix) && return view
    isempty(view) && return InsideView(view.indices, 1, 0, view.depth + length(prefix))

    first = view.first
    while first <= view.last && !_matches(view, view.indices[first], prefix)
        first += 1
    end
    last = first
    while last <= view.last && _matches(view, view.indices[last], prefix)
        last += 1
    end
    InsideView(view.indices, first, last - 1, view.depth + length(prefix))
end

child_view(view::InsideView, prefix::Tuple) = child_view(Val(true), view, prefix)
child_view(::Nothing, prefix::Tuple) = nothing
child_view(indices::AbstractVector, prefix::Tuple) = child_view(Val(true), InsideView(indices), prefix)

function child_view(::Val{true}, view::InsideView, first::Int, last::Int, depth::Int)
    InsideView(view.indices, first, last, depth)
end

child_view(view::InsideView, first::Int, last::Int, depth::Int) = child_view(Val(true), view, first, last, depth)

function has_other(view::InsideView, obstruction_indices::Tuple)
    for index in view
        index != obstruction_indices && return true
    end
    false
end

has_other(indices::AbstractVector, obstruction_indices::Tuple) = has_other(InsideView(indices), obstruction_indices)
has_other(::Nothing, ::Tuple) = false

end
