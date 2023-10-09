module RandomDistribution
import StaticArrays: SVector
import Distributions

function draw_radii(box_size::SVector{N, <:Real}, target_density::Real, distribution::Distributions.ContinuousUnivariateDistribution, allowed_range::Tuple{<:Real, <:Real}) where {N}
    if target_density > 1 || target_density < 0
        error("Target density should be a number between 0 and 1")
    end
    @assert allowed_range[2] > allowed_range[1]
    volume = 0.
    target_volume = prod(box_size) * target_density
    samples = Float64[]
    while volume < target_volume
        new_value = rand(distribution)
        if new_value < allowed_range[1] || new_value > allowed_range[2]
            continue 
        end
        push!(samples, new_value)
        if N == 3
            volume += 4//3 * π * new_value^3
        elseif N == 2
            volume += π * new_value^2
        elseif N == 1
            volume += 2 * new_value
        end
    end
    samples
end


function eliminate_overlap(positions::AbstractVector{SVector{N, Float64}}, radii::AbstractVector{Float64}, box_size::SVector{N, Float64}; max_iter=1000) where {N}
    @assert length(positions) == length(radii)
    moved = true
    half_box_size = box_size ./ 2

    step_size = 1.

    for _ in 1:max_iter
        moved = false
        for i in 1:length(positions)
            for j in 1:(i-1)
                displacement = map((dp, hbs) -> mod(dp + hbs, hbs * 2) - hbs, positions[i] - positions[j], half_box_size)
                distsq = sum(p->p*p, displacement)
                total_radius = radii[i] + radii[j]
                if distsq < (total_radius * total_radius)
                    moved = true
                    dist = sqrt(distsq)
                    total_movement = max((total_radius - dist) * (1 + step_size), 0.01 * total_radius)
                    # move larger cylinder less
                    positions[i] = positions[i] .+ displacement .* (total_movement * radii[j] / (total_radius * dist))
                    positions[j] = positions[j] .- displacement .* (total_movement * radii[i] / (total_radius * dist))
                end
            end
        end
        step_size *= 0.8
        if !moved
            return [map((p, hbs) -> mod(p + hbs, hbs * 2) - hbs, pos, half_box_size) for pos in positions]
        end
    end
    error("No solution without overlap found after max_iter=$max_iter iterations")
end


"""
    random_positions_radii(box_size, target_density, n_dimensions; distribution=Gamma, mean=1., variance=1., max_iter=1000, min_radius=0.1, max_radius=Inf)

Randomly distributes circles or spheres in space.

Arguments:
- `box_size`: Size of the infinitely repeating box of random positions
- `target_density`: Final density of the circles/spheres. This density will only be approximately reached
- `n_dimensions`: dimensionality of the space (2 for cicles; 3 for spheres)
- `distribution`: distribution from which the radii are drawn (from [Distributions.jl](https://juliastats.org/Distributions.jl/stable/))
- `mean`: mean of the gamma distribution (ignored if `distribution` explicitly set)
- `variance`: variance of the gamma distribution (ignored if `distribution` explicitly set)
- `max_iter`: maximum number of iterations to try to prevent the circles/spheres from overlapping. An error is raised if they still overlap after this number of iterations.
- `min_radius`: samples from the distribution below this radius will be rejected (in um).
- `max_radius`: samples from the distribution above this radius will be rejected (in um).
"""
function random_positions_radii(
    box_size, target_density::Real, ndim::Int;
    distribution=nothing, mean=1., variance=1., max_iter=1000,
    min_radius=0.1, max_radius=Inf,
)
    if iszero(variance)
        distribution = Distributions.Normal(mean, 0.)
    elseif isnothing(distribution)
        a = mean * mean / variance
        scale = variance/mean
        distribution = Distributions.Gamma(a, scale)
    end
    if isa(box_size, Real)
        box_size = SVector{ndim, Float64}(fill(box_size, ndim))
    else
        box_size = SVector{ndim, Float64}(box_size)
    end
    target_density = Float64(target_density)
    radii = draw_radii(box_size, target_density, distribution, (min_radius, max_radius))
    start_positions = [SVector{ndim, Float64}((rand(ndim) .- 0.5) .* box_size) for _ in 1:length(radii)]
    return eliminate_overlap(start_positions, radii, box_size; max_iter=max_iter), radii
end

end