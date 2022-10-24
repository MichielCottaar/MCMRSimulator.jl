@testset "Collision tests" begin
    function compare(ms1 :: AbstractVector{<:mr.Movement}, ms2 :: AbstractVector{<:mr.Movement})
        @test length(ms1) == length(ms2)
        for (m1, m2) in zip(ms1, ms2)
            @test m1.origin ≈ m2.origin atol=1e-7 rtol=1e-3
            @test m1.destination ≈ m2.destination atol=1e-7 rtol=1e-3
            @test m1.timestep ≈ m2.timestep atol=1e-7 rtol=1e-3
        end
    end
    @testset "Wall reflections" begin
        @testset "Hitting vertical wall directly" begin
            res = mr.correct_collisions(
                mr.Movement(SA[0, 0, 0], SA[3, 0, 0], 3.),
                mr.walls(rotation=:x, positions=[1])
            )
            @test length(res) == 2
            @test res[1].origin ≈ SA[0, 0, 0]
            @test res[1].destination ≈ SA[1, 0, 0]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA[1, 0, 0]
            @test res[2].destination ≈ SA[-1, 0, 0]
            @test res[2].timestep ≈ 2.
        end
        @testset "Hitting vertical wall under angle" begin
            res = mr.correct_collisions(
                mr.Movement(SA[0, 0, 0], SA[3, 6, 0], Float(3.)),
                mr.walls(rotation=:x, positions=[1])
            )
            @test length(res) == 2
            @test res[1].origin ≈ SA[0, 0, 0]
            @test res[1].destination ≈ SA[1, 2, 0]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA[1, 2, 0]
            @test res[2].destination ≈ SA[-1, 6, 0]
            @test res[2].timestep ≈ 2.
        end
        @testset "Missing vertical wall" begin
            res = mr.correct_collisions(
                mr.Movement(SA[0, 0, 0], SA[0, 2, 0], Float(4.)),
                mr.walls(rotation=:x, positions=[1])
            )
            @test length(res) == 1
            @test res[1].origin == SA[0, 0, 0]
            @test res[1].destination == SA[0, 2, 0]
            @test res[1].timestep == 4.
        end
        @testset "Hitting two vertical walls" begin
            res = mr.correct_collisions(
                mr.Movement(SA[0, 0, 0], SA[6, 0, 12], Float(30.)),
                mr.walls(rotation=:x, positions=[-1, 1])
            )
            @test length(res) == 4
            @test res[1].origin ≈ SA[0, 0, 0]
            @test res[1].destination ≈ SA[1, 0, 2]
            @test res[1].timestep ≈ 5.
            @test res[2].origin ≈ SA[1, 0, 2]
            @test res[2].destination ≈ SA[-1, 0, 6]
            @test res[2].timestep ≈ 10.
            @test res[3].origin ≈ SA[-1, 0, 6]
            @test res[3].destination ≈ SA[1, 0, 10]
            @test res[3].timestep ≈ 10.
            @test res[4].origin ≈ SA[1, 0, 10]
            @test res[4].destination ≈ SA[0, 0, 12]
            @test res[4].timestep ≈ 5.
        end
        @testset "Hitting vertical and horizontal walls" begin
            res = mr.correct_collisions(
                mr.Movement([-1, 0, 0], [2, 3, 3], 3.),
                [mr.walls(rotation=:x, positions=1.), mr.walls(rotation=:y, positions=1.)]
            )
            @test length(res) == 3
            @test res[1].origin ≈ SA[-1, 0, 0]
            @test res[1].destination ≈ SA[0, 1, 1]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA[0, 1, 1]
            @test res[2].destination ≈ SA[1, 0, 2]
            @test res[2].timestep ≈ 1.
            @test res[3].origin ≈ SA[1, 0, 2]
            @test res[3].destination ≈ SA[0, -1, 3]
            @test res[3].timestep ≈ 1.
        end
        @testset "Hitting a corner" begin
            res = mr.correct_collisions(
                mr.Movement(SA[0, 0, 0], SA[3, 3, 3], Float(3.)),
                [mr.walls(rotation=:x, positions=1.), mr.walls(rotation=:y, positions=1.)]
            )
            @test length(res) == 3
            @test res[1].origin ≈ SA[0, 0, 0]
            @test res[1].destination ≈ SA[1, 1, 1]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA[1, 1, 1]
            @test res[2].destination ≈ SA[1, 1, 1]
            @test res[2].timestep ≈ 0. atol=1e-5
            @test res[3].origin ≈ SA[1, 1, 1]
            @test res[3].destination ≈ SA[-1, -1, 3]
            @test res[3].timestep ≈ 2.
        end
        @testset "Hitting diagonal wall" begin
            res = mr.correct_collisions(
                mr.Movement([1, 0, 0], [1, 3, 0], 6.),
                mr.walls(rotation=[1, 1, 0], positions=sqrt(2)),
            )
            @test length(res) == 2
            @test res[1].origin ≈ SA[1, 0, 0]
            @test res[1].destination ≈ SA[1, 1, 0]
            @test res[1].timestep ≈ 2.
            @test res[2].origin ≈ SA[1, 1, 0]
            @test res[2].destination ≈ SA[-1, 1, 0]
            @test res[2].timestep ≈ 4.
        end
    end
    @testset "Sphere reflections" begin
        @testset "Remain within sphere" begin
            res = mr.correct_collisions(
                mr.Movement([0, 0, 1.5], [6, 8, 6], 6.),
                mr.spheres(1., positions=[0, 0, 2])
            )
            final = res[end].destination
            radius = norm(final .- SA[0, 0, 2])
            @assert radius <= 1.
        end
    end
    @testset "Cylinder reflections" begin
        @testset "Within cylinder along radial line" begin
            res = mr.correct_collisions(
                mr.Movement(SA[0, 0, 0], SA[6, 0, 6], 6.),
                [mr.cylinders(1., rotation=:z)]
            )
            @test length(res) == 4
            @test res[1].origin ≈ SA[0, 0, 0]
            @test res[1].destination ≈ SA[1, 0, 1]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA[1, 0, 1]
            @test res[2].destination ≈ SA[-1, 0, 3]
            @test res[2].timestep ≈ 2.
            @test res[3].origin ≈ SA[-1, 0, 3]
            @test res[3].destination ≈ SA[1, 0, 5]
            @test res[3].timestep ≈ 2.
            @test res[4].origin ≈ SA[1, 0, 5]
            @test res[4].destination ≈ SA[0, 0, 6]
            @test res[4].timestep ≈ 1.
        end
        @testset "90 degree bounces within vertical cylinder" begin
            res = mr.correct_collisions(
                mr.Movement(SA[0, 1, 0], SA[10, 1, 0], 10.),
                mr.cylinders(sqrt(2), rotation=:z, positions=[0, 0])
            )
            compare(res, [
                mr.Movement(SA[0, 1, 0], SA[1, 1, 0], 1.),
                mr.Movement(SA[1, 1, 0], SA[1, -1, 0], 2.),
                mr.Movement(SA[1, -1, 0], SA[-1, -1, 0], 2.),
                mr.Movement(SA[-1, -1, 0], SA[-1, 1, 0], 2.),
                mr.Movement(SA[-1, 1, 0], SA[1, 1, 0], 2.),
                mr.Movement(SA[1, 1, 0], SA[1, 0, 0], 1.),
            ])
        end
        @testset "90 degree bounces from outside vertical cylinder" begin
            res = mr.correct_collisions(
                mr.Movement(SA[-2, 1, 0], SA[2, 1, 0], 4.),
                [mr.cylinders(sqrt(2), rotation=:z, positions=[0, 0])]
            )
            compare(res, [
                mr.Movement(SA[-2, 1, 0], SA[-1, 1, 0], 1.),
                mr.Movement(SA[-1, 1, 0], SA[-1, 4, 0], 3.),
            ])
        end
        @testset "Remain within angled cylinder" begin
            orient = SA[1, 2, sqrt(3)]
            cylinder = mr.cylinders(2.3, rotation=orient)
            res = mr.correct_collisions(
                mr.Movement(SA[0, 0.5, 0.3], SA[-30, 50, 10], 40.),
                [cylinder]
            )
            final = res[end].destination
            radius = norm(final .- (orient ⋅ final) * orient / norm(orient) ^ 2)
            @test radius <= 2.3
            @test mr.isinside(cylinder, final)
        end
        @testset "Remain within distant cylinder" begin
            Random.seed!(1234)
            geometry = mr.cylinders([0.8, 0.9], repeats=[2., 2.])
            c2 = mr.Cylinder(0.8)
            spin = mr.Spin(position=SA[200., 200., 0.])
            @test mr.isinside(geometry, spin)
            inside = true
            for _ in 1:100
                spin = mr.draw_step(spin, Float(3.), Float(0.5), tuple(geometry))
                inside &= mr.isinside(geometry, spin)
            end
            @test inside
        end
    end
    @testset "Ray-grid intersections with undefined grid" begin
        function tcompare(t1, t2)
            @test length(t1) == length(t2)
            for (e1, e2) in zip(t1, t2)
                @test e1 ≈ e2
            end
        end
        res = collect(mr.ray_grid_intersections(SA[0.5, 0.5, 0.5], SA[0.5, 0.5, 3.5]))
        tcompare(res[1], ([0, 0, 0], 0., [0.5, 0.5, 0.5], 1/6, [0.5, 0.5, 1.]))
        tcompare(res[2], ([0, 0, 1], 1/6, [0.5, 0.5, 0.], 1/2, [0.5, 0.5, 1.]))
        tcompare(res[3], ([0, 0, 2], 1/2, [0.5, 0.5, 0.], 5/6, [0.5, 0.5, 1.]))
        tcompare(res[4], ([0, 0, 3], 5/6, [0.5, 0.5, 0.], 1., [0.5, 0.5, 0.5]))
    end
    @testset "Ray-grid intersections with defined grid" begin
        function tcompare(t1, t2)
            @test length(t1) == length(t2)
            for (e1, e2) in zip(t1, t2)
                @test e1 ≈ e2
            end
        end
        grid = mr.GridShape(mr.BoundingBox([0, 2, -Inf], [2, 12, Inf]), [4, 10, 1])
        @test mr.project(SA[0.25, 3.5, 0.5], grid) == SA[1.5, 2.5, 1.5]
        @test mr.project(SA[0.75, 5.5, 192.5], grid) == SA[2.5, 4.5, 1.5]
        res = collect(mr.ray_grid_intersections(grid, SA[0.25, 3.5, 0.5], SA[0.75, 5.5, 192.5]))
        tcompare(res[1], ([1, 2, 1], 0., [0.5, 0.5, 0.5], 1/4, [3/4, 1, 0.5]))
        tcompare(res[2], ([1, 3, 1], 1/4, [3/4, 0, 0.5], 1/2, [1, 0.5, 0.5]))
        tcompare(res[3], ([2, 3, 1], 1/2, [0, 0.5, 0.5], 3/4, [1/4, 1, 0.5]))
        tcompare(res[4], ([2, 4, 1], 3/4, [1/4, 0, 0.5], 1, [0.5, 0.5, 0.5]))
    end
    @testset "Ray-grid intersections with low-dimensional grid" begin
        function tcompare(t1, t2)
            @test length(t1) == length(t2)
            for (e1, e2) in zip(t1, t2)
                @test e1 ≈ e2
            end
        end
        grid = mr.GridShape(mr.BoundingBox([0, 2], [2, 12]), [4, 10])
        @test mr.project(SA[0.25, 3.5], grid) == SA[1.5, 2.5]
        @test mr.project(SA[0.75, 5.5], grid) == SA[2.5, 4.5]
        res = collect(mr.ray_grid_intersections(grid, SA[0.25, 3.5], SA[0.75, 5.5]))
        tcompare(res[1], ([1, 2], 0., [0.5, 0.5], 1/4, [3/4, 1]))
        tcompare(res[2], ([1, 3], 1/4, [3/4, 0], 1/2, [1, 0.5]))
        tcompare(res[3], ([2, 3], 1/2, [0, 0.5], 3/4, [1/4, 1]))
        tcompare(res[4], ([2, 4], 3/4, [1/4, 0], 1, [0.5, 0.5]))
    end
    @testset "Reflections on planes of cylinders" begin
        @testset "Bounce between four cylinders" begin
            cylinders = mr.cylinders(
                sqrt(2), repeats=[3, 4]
            )
            res = mr.correct_collisions(
                mr.Movement(SA[1, 2, 0], SA[1, 11, 9], 9.),
                [cylinders]
            )
            compare(res, [
                mr.Movement(SA[1, 2, 0], SA[1, 3, 1], 1.),
                mr.Movement(SA[1, 3, 1], SA[2, 3, 2], 1.),
                mr.Movement(SA[2, 3, 2], SA[2, 1, 4], 2.),
                mr.Movement(SA[2, 1, 4], SA[1, 1, 5], 1.),
                mr.Movement(SA[1, 1, 5], SA[1, 3, 7], 2.),
                mr.Movement(SA[1, 3, 7], SA[2, 3, 8], 1.),
                mr.Movement(SA[2, 3, 8], SA[2, 2, 9], 1.),
            ])
        end
        @testset "Travel through many repeats between bounces" begin
            cylinders = mr.cylinders(
                1, repeats=[2, 4]
            )
            res = mr.correct_collisions(
                mr.Movement([1, 2, 0], [13, 6, 4], 12.),
                [cylinders]
            )
            compare(res, [
                mr.Movement(SA[1, 2, 0], SA[4, 3, 1], 3.),
                mr.Movement(SA[4, 3, 1], SA[10, 1, 3], 6.),
                mr.Movement(SA[10, 1, 3], SA[13, 2, 4], 3.),
            ])
        end
        @testset "Test lots of particles still don't cross compartments" begin
            for geometry in (
                mr.cylinders([0.8, 0.9], positions=[[0, 0], [0, 0]], repeats=[2, 2]),
                mr.spheres([0.8, 0.9], positions=[[0, 0, 0], [2, 0, 2]], repeats=[2, 2, 2]),
            )
                sequence = mr.perfect_dwi(bval=2.)

                snap = mr.Snapshot(300);

                simulation = mr.Simulation([sequence], diffusivity=3., geometry=geometry);

                before = mr.isinside(geometry, snap);
                after = mr.isinside(geometry, mr.evolve(snap, simulation, 200.))
                switched = sum(xor.(before, after))
                @test switched == 0
            end
        end
        @testset "Test that we remain between two walls" begin
            Random.seed!(1234)
            walls = mr.walls(positions=[0, 1])
            snap = mr.Snapshot([mr.Spin(position=rand(3)) for _ in 1:3000])
            sequence = mr.perfect_dwi(bval=2.)
            simulation = mr.Simulation([sequence]; geometry=walls, diffusivity=3.);

            final = mr.evolve(snap, simulation, 200.)
            xfinal = [s.position[1] for s in final.spins]
            @test all(xfinal .>= 0.)
            @test all(xfinal .<= 1.)
        end
    end
    @testset "Test spiral collision detection" begin
        spiral = mr.spirals(1., 10., thickness=1., inner_cylinder=false)

        function test(origin, destination, distance, normal)
            m = mr.Movement(origin, destination, 1)
            c = mr.detect_collision(m, spiral, mr.empty_collision)
            @test c.distance ≈ distance
            @test (c.normal ./ norm(c.normal)) ≈ (mr.PosVector(normal) ./ norm(normal))
        end

        test([0.5, 0.5, 0], [1.5, 1.5, 0], (1.125 / sqrt(2)) - 0.5, [-1, -1, 0])
        test([0., 0., 0], [1., 1., 0], 1.125 / sqrt(2), [-1, -1, 0])
        test([0, 1, 0], [0, 2, 0], 0.25, [0, -1, 0])
        test([-1, 0, 0], [-2, 0, 0], 0.5, [1, 0, 0])
        test([1, 1, 0], [0.5, 0.5, 0], (1 - 1.125 / sqrt(2)) * 2, [1, 1, 0])
        test([1, 1, 0], [0., 0., 0], 1 - 1.125 / sqrt(2), [1, 1, 0])
        test([0, 2, 0], [0, 1, 0], 0.75, [0, 1, 0])
        test([-2, 0, 0], [-1, 0, 0], 0.5, [-1, 0, 0])

        test([0, -1, 0], [4, 1, 0], 0.5, [-1, 0, 0])
        m = mr.Movement([0, 1, 0], [4, -1, 0], 1)
        c = mr.detect_collision(m, spiral, mr.empty_collision)
        @test 0 < c.distance < 0.25

        test([11, 1, 0], [7, -1, 0], 0.5, [1, 0, 0])
        m = mr.Movement([11, -1, 0], [7, 1, 0], 1)
        c = mr.detect_collision(m, spiral, mr.empty_collision)
        @test 0.25 < c.distance < 0.3
    end
    @testset "Spiral leakage detection" begin
        spiral = mr.spirals(0.8, 1., inner_cylinder=true, outer_cylinder=true)
        spin = mr.Spin(position=[0.9, 0., 0.])
        theta = mr.spiral_theta(spiral.obstructions[1], SVector{2}(spin.position[1:2]))
        @show spin.position
        @show theta
        for _ in 1:100000
            prev_theta = theta
            spin = mr.draw_step(spin, 1., 0.01, [spiral])
            @show spin.position
            theta = mr.spiral_theta(spiral.obstructions[1], SVector{2}(spin.position[1:2]))
            @show theta
            if abs(theta - prev_theta) > π
                error("Leaked through spiral!")
            end
            if mr.isinside(spiral, spin) != 1
                error("Leaked through cylinders!")
            end
        end
    end
end