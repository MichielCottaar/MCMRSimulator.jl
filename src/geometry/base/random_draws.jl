function draw_radii(box_size::SVector{N, <:Real}, target_density::Real, distribution::Distributions.ContinuousUnivariateDistribution) where {N}
    if target_density > 1 || target_density < 0
        error("Target density should be a number between 0 and 1")
    end
    volume = 0.
    target_volume = prod(box_size) * target_density
    samples = Float[]
    while volume < target_volume
        new_value = rand(distribution)
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


function eliminate_overlap(positions::AbstractVector{SVector{N, Float}}, radii::AbstractVector{Float}, box_size::SVector{N, Float}; max_iter=1000) where {N}
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
                    total_movement = max((total_radius - sqrt(distsq)) * (1 + step_size), 0.01 * total_radius)
                    #println(i, " ", j, " ", total_movement)
                    # move larger cylinder less
                    positions[i] = positions[i] .+ displacement .* (total_movement * radii[j] / total_radius)
                    positions[j] = positions[j] .- displacement .* (total_movement * radii[i] / total_radius)
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


function random_positions_radii(
    box_size, target_density::Real, ndim::Int;
    distribution=nothing, mean=1., variance=1., max_iter=1000
)
    if isnothing(distribution)
        a = mean * mean / variance
        scale = variance/mean
        distribution = Distributions.Gamma(a, scale)
    end
    if isa(box_size, Real)
        box_size = SVector{ndim, Float}(fill(box_size, ndim))
    else
        box_size = SVector{ndim, Float}(box_size)
    end
    target_density = Float(target_density)
    radii = draw_radii(box_size, target_density, distribution)
    start_positions = [SVector{ndim, Float}((rand(ndim) .- 0.5) .* box_size) for _ in 1:length(radii)]
    return eliminate_overlap(start_positions, radii, box_size; max_iter=max_iter), radii
end