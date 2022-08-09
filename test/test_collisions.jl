@testset "Collision tests" begin
    function compare(ms1 :: AbstractVector{<:mr.Movement}, ms2 :: AbstractVector{<:mr.Movement})
        @test length(ms1) == length(ms2)
        for (m1, m2) in zip(ms1, ms2)
            @test m1.origin ≈ m2.origin atol=1e-12 rtol=1e-6
            @test m1.destination ≈ m2.destination atol=1e-12 rtol=1e-6
            @test m1.timestep ≈ m2.timestep atol=1e-12 rtol=1e-6
        end
    end
    @testset "Wall reflections" begin
        @testset "Hitting vertical wall directly" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 0, 0], SA_F64[3, 0, 0], 3.),
                mr.Wall(:x, 1.),
            )
            @test length(res) == 2
            @test res[1].origin ≈ SA_F64[0, 0, 0]
            @test res[1].destination ≈ SA_F64[1, 0, 0]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA_F64[1, 0, 0]
            @test res[2].destination ≈ SA_F64[-1, 0, 0]
            @test res[2].timestep ≈ 2.
        end
        @testset "Hitting vertical wall under angle" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 0, 0], SA_F64[3, 6, 0], 3.),
                mr.Wall(:x, 1.),
            )
            @test length(res) == 2
            @test res[1].origin ≈ SA_F64[0, 0, 0]
            @test res[1].destination ≈ SA_F64[1, 2, 0]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA_F64[1, 2, 0]
            @test res[2].destination ≈ SA_F64[-1, 6, 0]
            @test res[2].timestep ≈ 2.
        end
        @testset "Missing vertical wall" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 0, 0], SA_F64[0, 2, 0], 4.),
                mr.Wall(:x, 1.),
            )
            @test length(res) == 1
            @test res[1].origin == SA_F64[0, 0, 0]
            @test res[1].destination == SA_F64[0, 2, 0]
            @test res[1].timestep == 4.
        end
        @testset "Hitting two vertical walls" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 0, 0], SA_F64[6, 0, 12], 30.),
                [mr.Wall(:x, 1.), mr.Wall(:x, -1.)],
            )
            @test length(res) == 4
            @test res[1].origin ≈ SA_F64[0, 0, 0]
            @test res[1].destination ≈ SA_F64[1, 0, 2]
            @test res[1].timestep ≈ 5.
            @test res[2].origin ≈ SA_F64[1, 0, 2]
            @test res[2].destination ≈ SA_F64[-1, 0, 6]
            @test res[2].timestep ≈ 10.
            @test res[3].origin ≈ SA_F64[-1, 0, 6]
            @test res[3].destination ≈ SA_F64[1, 0, 10]
            @test res[3].timestep ≈ 10.
            @test res[4].origin ≈ SA_F64[1, 0, 10]
            @test res[4].destination ≈ SA_F64[0, 0, 12]
            @test res[4].timestep ≈ 5.
        end
        @testset "Hitting vertical and horizontal walls" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[-1, 0, 0], SA_F64[2, 3, 3], 3.),
                [mr.Wall(:x, 1.), mr.Wall(:y, 1.)],
            )
            @test length(res) == 3
            @test res[1].origin ≈ SA_F64[-1, 0, 0]
            @test res[1].destination ≈ SA_F64[0, 1, 1]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA_F64[0, 1, 1]
            @test res[2].destination ≈ SA_F64[1, 0, 2]
            @test res[2].timestep ≈ 1.
            @test res[3].origin ≈ SA_F64[1, 0, 2]
            @test res[3].destination ≈ SA_F64[0, -1, 3]
            @test res[3].timestep ≈ 1.
        end
        @testset "Hitting a corner" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 0, 0], SA_F64[3, 3, 3], 3.),
                [mr.Wall(:x, 1.), mr.Wall(:y, 1.)],
            )
            @test length(res) == 3
            @test res[1].origin ≈ SA_F64[0, 0, 0]
            @test res[1].destination ≈ SA_F64[1, 1, 1]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA_F64[1, 1, 1]
            @test res[2].destination ≈ SA_F64[1, 1, 1]
            @test res[2].timestep ≈ 0. atol=1e-8
            @test res[3].origin ≈ SA_F64[1, 1, 1]
            @test res[3].destination ≈ SA_F64[-1, -1, 3]
            @test res[3].timestep ≈ 2.
        end
        @testset "Hitting diagonal wall" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[1, 0, 0], SA_F64[1, 3, 0], 6.),
                mr.Wall(SA_F64[1, 1, 0], 1.),
            )
            @test length(res) == 2
            @test res[1].origin ≈ SA_F64[1, 0, 0]
            @test res[1].destination ≈ SA_F64[1, 1, 0]
            @test res[1].timestep ≈ 2.
            @test res[2].origin ≈ SA_F64[1, 1, 0]
            @test res[2].destination ≈ SA_F64[-1, 1, 0]
            @test res[2].timestep ≈ 4.
        end
    end
    @testset "Sphere reflections" begin
        @testset "Remain within sphere" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 0, 1.5], SA_F64[6, 8, 6], 6.),
                [mr.Sphere(1., SA_F64[0, 0, 2])]
            )
            final = res[end].destination
            radius = norm(final .- SA_F64[0, 0, 2])
            @assert radius <= 1.
        end
    end
    @testset "Cylinder reflections" begin
        @testset "Within cylinder along radial line" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 0, 0], SA_F64[6, 0, 6], 6.),
                [mr.Cylinder(1., :z, SA_F64[0, 0, 2])]
            )
            @test length(res) == 4
            @test res[1].origin ≈ SA_F64[0, 0, 0]
            @test res[1].destination ≈ SA_F64[1, 0, 1]
            @test res[1].timestep ≈ 1.
            @test res[2].origin ≈ SA_F64[1, 0, 1]
            @test res[2].destination ≈ SA_F64[-1, 0, 3]
            @test res[2].timestep ≈ 2.
            @test res[3].origin ≈ SA_F64[-1, 0, 3]
            @test res[3].destination ≈ SA_F64[1, 0, 5]
            @test res[3].timestep ≈ 2.
            @test res[4].origin ≈ SA_F64[1, 0, 5]
            @test res[4].destination ≈ SA_F64[0, 0, 6]
            @test res[4].timestep ≈ 1.
        end
        @testset "90 degree bounces within vertical cylinder" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 1, 0], SA_F64[10, 1, 0], 10.),
                [mr.Cylinder(sqrt(2), :z, SA_F64[0, 0, 2])]
            )
            compare(res, [
                mr.Movement(SA_F64[0, 1, 0], SA_F64[1, 1, 0], 1.),
                mr.Movement(SA_F64[1, 1, 0], SA_F64[1, -1, 0], 2.),
                mr.Movement(SA_F64[1, -1, 0], SA_F64[-1, -1, 0], 2.),
                mr.Movement(SA_F64[-1, -1, 0], SA_F64[-1, 1, 0], 2.),
                mr.Movement(SA_F64[-1, 1, 0], SA_F64[1, 1, 0], 2.),
                mr.Movement(SA_F64[1, 1, 0], SA_F64[1, 0, 0], 1.),
            ])
        end
        @testset "90 degree bounces from outside vertical cylinder" begin
            res = mr.correct_collisions(
                mr.Movement(SA_F64[-2, 1, 0], SA_F64[2, 1, 0], 4.),
                [mr.Cylinder(sqrt(2), :z, SA_F64[0, 0, 2])]
            )
            compare(res, [
                mr.Movement(SA_F64[-2, 1, 0], SA_F64[-1, 1, 0], 1.),
                mr.Movement(SA_F64[-1, 1, 0], SA_F64[-1, 4, 0], 3.),
            ])
        end
        @testset "Remain within angled cylinder" begin
            orient = SA_F64[1, 2, sqrt(3)]
            cylinder = mr.Cylinder(2.3, orient, SA_F64[0, 0, 0])
            res = mr.correct_collisions(
                mr.Movement(SA_F64[0, 0.5, 0.3], SA_F64[-30, 50, 10], 40.),
                [cylinder]
            )
            final = res[end].destination
            radius = norm(final .- (orient ⋅ final) * orient / norm(orient) ^ 2)
            @test radius <= 2.3
            @test mr.isinside(final, cylinder)
        end
        @testset "Remain within distant cylinder" begin
            Random.seed!(1234)
            cylinder = mr.Cylinder(0.9)
            c2 = mr.Cylinder(0.8)
            geometry = SVector{1}(mr.Repeated([cylinder, c2], [2., 2., Inf]))
            position = SA_F64[200., 200., 0.]
            @test mr.isinside(position, geometry)
            for _ in 1:100
                position = mr.draw_step(position, 3., 0.5, geometry)
                @test mr.isinside(position, geometry)
            end
        end
    end
    @testset "Ray-grid intersections" begin
        function tcompare(t1, t2)
            @test length(t1) == length(t2)
            for (e1, e2) in zip(t1, t2)
                @test e1 ≈ e2
            end
        end
        res = collect(mr.ray_grid_intersections(SA_F64[0.5, 0.5, 0.5], SA_F64[0.5, 0.5, 3.5]))
        tcompare(res[1], ([0, 0, 0], 0., [0.5, 0.5, 0.5], 1/6, [0.5, 0.5, 1.]))
        tcompare(res[2], ([0, 0, 1], 1/6, [0.5, 0.5, 0.], 1/2, [0.5, 0.5, 1.]))
        tcompare(res[3], ([0, 0, 2], 1/2, [0.5, 0.5, 0.], 5/6, [0.5, 0.5, 1.]))
        tcompare(res[4], ([0, 0, 3], 5/6, [0.5, 0.5, 0.], 1., [0.5, 0.5, 0.5]))
    end
    @testset "Reflections on planes of cylinders" begin
        @testset "Bounce between four cylinders" begin
            cylinders = mr.cylinder_plane(
                sqrt(2), repeatx=3, repeaty=4
            )
            res = mr.correct_collisions(
                mr.Movement(SA_F64[1, 2, 0], SA_F64[1, 11, 9], 9.),
                [cylinders]
            )
            compare(res, [
                mr.Movement(SA_F64[1, 2, 0], SA_F64[1, 3, 1], 1.),
                mr.Movement(SA_F64[1, 3, 1], SA_F64[2, 3, 2], 1.),
                mr.Movement(SA_F64[2, 3, 2], SA_F64[2, 1, 4], 2.),
                mr.Movement(SA_F64[2, 1, 4], SA_F64[1, 1, 5], 1.),
                mr.Movement(SA_F64[1, 1, 5], SA_F64[1, 3, 7], 2.),
                mr.Movement(SA_F64[1, 3, 7], SA_F64[2, 3, 8], 1.),
                mr.Movement(SA_F64[2, 3, 8], SA_F64[2, 2, 9], 1.),
            ])
        end
        @testset "Travel through many repeats between bounces" begin
            cylinders = mr.cylinder_plane(
                1, repeatx=2, repeaty=4
            )
            res = mr.correct_collisions(
                mr.Movement(SA_F64[1, 2, 0], SA_F64[13, 6, 4], 12.),
                [cylinders]
            )
            compare(res, [
                mr.Movement(SA_F64[1, 2, 0], SA_F64[4, 3, 1], 3.),
                mr.Movement(SA_F64[4, 3, 1], SA_F64[10, 1, 3], 6.),
                mr.Movement(SA_F64[10, 1, 3], SA_F64[13, 2, 4], 3.),
            ])
        end
    end
end
