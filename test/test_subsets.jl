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
@testset "Test multiple geometries" begin
    geometry = [mr.Spheres(position=[0., 0., 0.], radius=1.), mr.Walls(position=[0., 1.2], density=1., dwell_time=1.)]
    sim = mr.Simulation([], geometry=geometry)
    N = 100000
    snapshot = mr.Snapshot(N, sim, 1.5)

    @test length(mr.get_subset(snapshot, sim)) == N
    @test length(mr.get_subset(snapshot, sim, bound=true)) ≈ 0.4 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=false)) ≈ 0.6 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=true, geometry_index=1)) == 0
    @test length(mr.get_subset(snapshot, sim, bound=false, geometry_index=1)) == N
    @test length(mr.get_subset(snapshot, sim, bound=true, geometry_index=2)) ≈ 0.4 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=false, geometry_index=2)) ≈ 0.6 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=true, geometry_index=2, obstruction_index=1)) ≈ 0.2 * N rtol = 0.1
    @test length(mr.get_subset(snapshot, sim, bound=true, geometry_index=2, obstruction_index=2)) ≈ 0.2 * N rtol = 0.1
    @test length(mr.get_subset(mr.get_subset(snapshot, sim, bound=true, geometry_index=2, obstruction_index=2), sim, inside=true)) == 0
    @test length(mr.get_subset(mr.get_subset(snapshot, sim, bound=true, geometry_index=2, obstruction_index=2), sim, inside=false)) ≈ 0.2 * N rtol=0.1

end

end