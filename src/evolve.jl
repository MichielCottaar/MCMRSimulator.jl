"""
    draw_step!(spin, sequence_parts, diffusivity, timestep, default_properties, geometry])

Updates the spin based on a random movement through the given geometry for a given `timestep`:
- draws the next location of the particle after `timestep` with given `diffusivity`.  
  This displacement will take into account the obstructions in `geometry`.
- The spin orientation will be affected by relaxation (see [`relax!`](@ref)) and potentially by magnetisation transfer during collisions.
"""
function draw_step!(spin :: Spin{N}, parts::SVector{N, SequencePart}, diffusivity :: Float, timestep :: Float, default_properties::GlobalProperties, geometry :: Geometry{0}) where {N}
    if iszero(timestep)
        return
    end
    relax!(spin, parts, geometry, 0, 1//2, default_properties.mri)
    @spin_rng spin begin
        spin.position += randn(PosVector) .* sqrt(2 * diffusivity * timestep)
    end
    relax!(spin, parts, geometry, 1//2, 1, default_properties.mri)
end

function draw_step!(spin :: Spin{N}, parts::SVector{N, SequencePart}, diffusivity :: Float, timestep :: Float, default_properties::GlobalProperties, geometry::Geometry) where {N}
    if iszero(timestep)
        return
    end
    current_pos = spin.position
    reflection = spin.stuck_to
    is_stuck = stuck(spin)
    current_time = zero(Float)

    found_solution = false
    @spin_rng spin begin
        if !stuck(spin)
            displacement = randn(PosVector) .* sqrt(2 * diffusivity * timestep)
            new_pos = current_pos + displacement
        end
        for _ in 1:10000
            if is_stuck
                td = dwell_time(stuck_to(spin), default_properties)
                time_stuck = -log(rand()) * td / timestep
                relax!(spin, parts, geometry, current_time, min(one(Float), current_time + time_stuck), default_properties.mri)
                current_time += time_stuck
                if current_time >= 1
                    found_solution = true
                    break
                end
                if iszero(Reflection(spin).timestep)
                    # special case where spin was generated as stuck on the surface
                    displacement = randn(PosVector) .* sqrt(2 * diffusivity * timestep * (1 - current_time))
                    if displacement ⋅ Collision(spin).normal < 0
                        displacement = -displacement
                    end
                else
                    displacement = direction(Reflection(spin), (1 - current_time) * timestep)
                end
                new_pos = current_pos + displacement

                reflection = spin.stuck_to
                spin.stuck_to = empty_reflection
                is_stuck = false
            end

            movement = Movement(current_pos, new_pos, one(Float))
            collision = detect_collision(
                movement,
                geometry,
                reflection.collision
            )

            use_distance = collision === empty_collision ? 1 : collision.distance
            next_time = current_time + (1 - current_time) * use_distance
            relax_pos_dist = rand() * use_distance
            spin.position = relax_pos_dist .* new_pos .+ (1 - relax_pos_dist) .* current_pos
            relax!(spin, parts, geometry, current_time, next_time, default_properties.mri)

            if collision === empty_collision
                spin.position = new_pos
                found_solution = true
                break
            end
            transfer!.(spin.orientations, correct_for_timestep(MT_fraction(collision.properties, default_properties), timestep))

            permeability_prob = correct_for_timestep(permeability(collision.properties, default_properties), timestep)
            passes_through = isone(permeability_prob) || !(iszero(permeability_prob) || rand() > permeability_prob)
            reflection = Reflection(collision, movement, timestep, next_time, passes_through)
            current_pos = spin.position = @. (movement.origin * (1 - collision.distance) + movement.destination * collision.distance)

            sd = surface_density(collision.properties, default_properties)
            if !iszero(sd) && rand() < stick_probability(sd, dwell_time(collision.properties, default_properties), diffusivity, timestep)
                spin.stuck_to = reflection
                is_stuck = true
            else
                new_pos = current_pos .+ direction(reflection)
            end

            current_time = next_time
        end
    end
    if !found_solution
        error("Bounced single particle for 10000 times in single step; terminating!")
    end
end


"""
    evolve_to_time(spin, simulation, current_time, new_time)

Evolve a single spin to the next time of interest.
This takes into account both random diffusion of the spin's position
and relaxation of the MR spin orientation.
It is used internally when evolving [`Simulation`](@ref) objects.
"""
function evolve_to_time!(
    spin::Spin{N}, simulation::Simulation{N}, parts::SVector{N, SequencePart}, current_time::Float, new_time::Float
) where {N}
    if current_time > new_time
        throw(DomainError("Spins cannot travel backwards in time"))
    end
    if new_time == current_time
        return spin
    end

    draw_step!(spin, parts, simulation.diffusivity, new_time - current_time, simulation.properties, simulation.geometry)
end

"""
    evolve_to_time(snapshot, simulation, new_time)

Evolves the full [`Snapshot`](@ref) through the [`Simulation`](@ref) to the given `new_time`.
Multi-threading is used to evolve multiple spins in parallel.
This is used internally when calling any of the snapshot evolution methods (e.g., [`evolve`](@ref)).
"""
function evolve_to_time(snapshot::Snapshot{N}, simulation::Simulation{N}, new_time::Float) where {N}
    current_time::Float = snapshot.time
    if new_time < current_time
        error("New requested time ($(new_time)) is less than current time ($(snapshot.time)). Simulator does not work backwards in time.")
    end
    spins::Vector{Spin{N}} = deepcopy.(snapshot.spins)

    times = propose_times(simulation, snapshot.time, new_time)

    # define next stopping times due to sequence, readout or times
    sequence_instants = Union{Nothing, InstantComponent}[next_instant(seq, current_time) for seq in simulation.sequences]
    sequence_times = MVector{N, Float}([isnothing(i) ? Inf : get_time(i) for i in sequence_instants])

    for next_time in times
        # evolve all spins to next interesting time
        parts = SequencePart.(simulation.sequences, current_time, next_time)
        Threads.@threads for spin in spins
            draw_step!(spin, parts, simulation.diffusivity, next_time - current_time, simulation.properties, simulation.geometry)
        end
        current_time = next_time

        # return final snapshot state
        if current_time != new_time && any(t -> t == current_time, sequence_times)
            components = SVector{N, Union{Nothing, InstantComponent}}([
                time == current_time ? instant : nothing 
                for (seq, instant, time) in zip(simulation.sequences, sequence_instants, sequence_times)
            ])
            apply!(components, spins)
            for (idx, ctime) in enumerate(sequence_times)
                if ctime == current_time
                    sequence = simulation.sequences[idx]
                    sequence_instants[idx] = next_instant(sequence, nextfloat(ctime))
                    sequence_times[idx] = get_time(sequence_instants[idx])
                end
            end
        end
    end
    return Snapshot(spins, current_time)
end