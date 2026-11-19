@testset "test_subsets.jl" begin

@testset "In empty geometry all spins are free and outside" begin
    for geometry in [[], mr.Walls(position=0., repeats=1.)]
        sim = mr.Simulation([], geometry=geometry)
        snapshot = mr.Snapshot(1000, sim)
        @test length(mr.get_subset(snapshot, sim)) == 1000
        @test length(mr.get_subset(snapshot, sim, bound=true)) == 0
        @test length(mr.get_subset(snapshot, sim, inside=true)) == 0
        @test length(mr.get_subset(snapshot, sim, bound=false)) == 1000
        @test length(mr.get_subset(snapshot, sim, inside=false)) == 1000
    end
end
@testset "Non-repeating geometry defines the default sampling box" begin
    geometry = mr.Spheres(position=[10., 20., 30.], radius=1.)
    simulation = mr.Simulation([], geometry=geometry)
    box = mr.BoundingBox(simulation.geometry)
    @test box.lower == [9., 19., 29.]
    @test box.upper == [11., 21., 31.]

    snapshot = mr.Snapshot(1000, simulation)
    @test all(all(mr.position(spin) .>= box.lower) && all(mr.position(spin) .<= box.upper) for spin in snapshot)
end
@testset "In repeated walled environment with surface density some spins are bound" begin
    geometry = mr.Walls(position=0., repeats=1., density=1., dwell_time=1.)
    sim = mr.Simulation([], geometry=geometry)
    snapshot = mr.Snapshot(100000, sim)
    @test length(mr.get_subset(snapshot, sim)) == 100000
    @test length(mr.get_subset(snapshot, sim, bound=true)) ≈ 50000 rtol=0.1
    @test length(mr.get_subset(snapshot, sim, inside=true)) == 0
    @test length(mr.get_subset(snapshot, sim, bound=false)) ≈ 50000 rtol=0.1
    @test length(mr.get_subset(snapshot, sim, inside=false)) == 100000
    bound_particles = mr.get_subset(snapshot, sim, bound=true)
    @test all(spin->iszero(mod(spin.position[1], 1.)), bound_particles)
    @test mr.stuck(sim, bound_particles[1])
    @test mr.stuck(sim, geometry, bound_particles[1])
    @test mr.stuck(sim, geometry[1], bound_particles[1])
    @test !mr.stuck(sim, mr.Walls(position=0.5, repeats=1.), bound_particles[1])
    @test all(mr.stuck(sim, geometry, mr.Snapshot(bound_particles)))
end
@testset "Testing inside of single sphere" begin
    geometry = mr.Spheres(position=[0., 0., 0.], radius=1.)
    sim = mr.Simulation([], geometry=geometry)
    snapshot = mr.Snapshot(100000, sim, 1.5)
    fraction = (4/3 * π) / 3^3
    @test length(mr.get_subset(snapshot, sim)) == 100000
    @test length(mr.get_subset(snapshot, sim, bound=true)) == 0
    @test length(mr.get_subset(snapshot, sim, inside=true)) ≈ fraction * 100000 rtol=0.1
    @test length(mr.get_subset(snapshot, sim, inside=true, bound=false)) ≈ fraction * 100000 rtol=0.1
    @test length(mr.get_subset(snapshot, sim, inside=true, bound=true)) == 0
    @test length(mr.get_subset(snapshot, sim, bound=false)) == 100000
    @test length(mr.get_subset(snapshot, sim, inside=false)) ≈ (1 - fraction) * 100000 rtol=0.1
    @test length(mr.get_subset(snapshot, sim, inside=false, bound=false)) ≈ (1 - fraction) * 100000 rtol=0.1
    @test length(mr.get_subset(snapshot, sim, inside=false, bound=true)) == 0
    inside_particles = mr.get_subset(snapshot, sim, inside=true)
    @test all(spin->norm(spin.position) < 1, inside_particles)
    outside_particles = mr.get_subset(snapshot, sim, inside=false)
    @test all(spin->norm(spin.position) > 1, outside_particles)
end

@testset "Use cached inside state for liminal geometry" begin
    Random.seed!(1234)
    geometry = mr.LiminalGeometry(
        geometries=[(1., mr.Spheres(radius=1.))],
        extracellular_fraction=0.5,
    )
    simulation = mr.Simulation([], geometry=geometry)
    snapshot = mr.Snapshot(100, simulation, 2.)
    expected_inside = count(!isempty, getfield.(snapshot.spins, :isinside))

    @test expected_inside > 0
    @test expected_inside < length(snapshot)
    @test length(mr.get_subset(snapshot, simulation, inside=true)) == expected_inside
    @test length(mr.get_subset(snapshot, simulation, inside=false)) ==
        length(snapshot) - expected_inside
end

@testset "Testing an explicit user geometry" begin
    sphere = mr.Spheres(position=[0.0, 0.0, 0.0], radius=1.0)
    simulation = mr.Simulation([], geometry=sphere)
    snapshot = mr.Snapshot([
        mr.Spin(position=[0.0, 0.0, 0.0]),
        mr.Spin(position=[1.0, 0.0, 0.0]),
        mr.Spin(position=[2.0, 0.0, 0.0]),
    ])

    @test length(mr.get_subset(snapshot, simulation, geometry=sphere, inside=true)) == 1
    @test length(mr.get_subset(snapshot, simulation, geometry=sphere, inside=false)) == 2
    @test length(mr.get_subset(snapshot, simulation, geometry=sphere, bound=true)) == 0
    @test length(mr.get_subset(snapshot, simulation, geometry=sphere, bound=false)) == 3

    duplicate_sphere = mr.Spheres(position=[0.0, 0.0, 0.0], radius=1.0)
    @test length(mr.get_subset(snapshot, simulation, geometry=duplicate_sphere, inside=true)) == 0
    @test length(mr.get_subset(snapshot, simulation, geometry=duplicate_sphere, inside=false)) == 3
end

@testset "Subset geometry can be an obstruction-group view" begin
    spheres = mr.Spheres(
        position=[[0., 0., 0.], [10., 0., 0.], [20., 0., 0.]],
        radius=1.,
    )
    simulation = mr.Simulation([], geometry=spheres)
    snapshot = mr.Snapshot([
        mr.Spin(position=[0., 0., 0.]),
        mr.Spin(position=[10., 0., 0.]),
        mr.Spin(position=[20., 0., 0.]),
    ])

    selected = mr.get_subset(snapshot, simulation; geometry=spheres[1:2], inside=true)
    @test length(selected) == 2
    @test all(spin -> spin.position[1] in (0., 10.), selected)
end

@testset "isinside matches simulation obstruction identities" begin
    geometry = mr.Spheres(
        position=[[0., 0., 0.], [10., 0., 0.], [20., 0., 0.]],
        radius=1.,
    )
    simulation = mr.Simulation([], geometry=geometry)
    snapshot = mr.Snapshot([
        mr.Spin(position=[0., 0., 0.]),
        mr.Spin(position=[10., 0., 0.]),
        mr.Spin(position=[20., 0., 0.]),
    ])

    @test mr.isinside(simulation, snapshot) == [1, 1, 1]
    @test mr.isinside(simulation, geometry[1:2], snapshot.spins[1]) == 1
    @test mr.isinside(simulation, geometry[1:2], snapshot.spins[3]) == 0
    @test mr.isinside(simulation, geometry[3:3], snapshot.spins[3]) == 1

    cached_snapshot = mr.Simulations._constrain_snapshot(snapshot, simulation)
    @test all(!isnothing(spin.property_indices) for spin in cached_snapshot)
    @test mr.isinside(simulation, cached_snapshot.spins[2]) == 1
    @test mr.isinside(simulation, geometry[1:2], cached_snapshot.spins[3]) == 0
end

@testset "Test multiple geometries" begin
    geometry = [mr.Spheres(position=[0., 0., 0.], radius=1.), mr.Walls(position=[0., 1.2], density=1., dwell_time=1.)]
    sim = mr.Simulation([], geometry=geometry)
    N = 100000
    snapshot = mr.Snapshot(N, sim, 1.5)

    @test length(mr.get_subset(snapshot, sim)) == N
    @test length(mr.get_subset(snapshot, sim, bound=true)) ≈ 0.4 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=false)) ≈ 0.6 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=true, geometry=geometry[1])) == 0
    @test length(mr.get_subset(snapshot, sim, bound=false, geometry=geometry[1])) == N
    @test length(mr.get_subset(snapshot, sim, bound=true, geometry=geometry[2])) ≈ 0.4 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=false, geometry=geometry[2])) ≈ 0.6 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=true, geometry=geometry[2][1])) ≈ 0.2 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=true, geometry=geometry[2][2])) ≈ 0.2 * N rtol = 0.1
    selected_snapshot = mr.get_subset(snapshot, sim, bound=true, geometry=geometry[2][2])
    @test length(mr.get_subset(selected_snapshot, sim, inside=true)) == 0
    @test length(mr.get_subset(selected_snapshot, sim, inside=false)) ≈ 0.2 * N rtol=0.1

end

end
