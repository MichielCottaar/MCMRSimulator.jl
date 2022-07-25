struct Movement
    origin :: SVector{3, Real}
    destination :: SVector{3, Real}
    timestep :: Real
end

draw_step(diffusivity :: Real, timestep :: Real) = sqrt(timestep * diffusivity) * @SVector randn(3)
draw_step(current_pos :: SVector, diffusivity :: Real, timestep :: Real) = Movement(
    current_pos,
    current_pos + draw_step(diffusivity, timestep),
    timestep
)


