@testset "test_collisions.jl" begin

    function compare(ms1 :: AbstractVector, ms2 :: AbstractVector)
        @test length(ms1) == length(ms2)
        for (m1, m2) in zip(ms1, ms2)
            @test m1 ≈ m2
        end
    end
    @testset "Wall reflections" begin
        @testset "Hitting vertical wall directly" begin
            res = correct_collisions(
                [0, 0, 0], [3, 0, 0],
                mr.Walls(rotation=:x, position=[1])
            )
            compare(res, [[0, 0, 0], [1, 0, 0], [-1, 0, 0]])
        end
        @testset "Hitting vertical wall under angle" begin
            res = correct_collisions(
                [0, 0, 0], [3, 6, 0],
                mr.Walls(rotation=:x, position=[1])
            )
            compare(res, [
                [0, 0, 0], [1, 2, 0], [-1, 6, 0]
            ])
        end
        @testset "Missing vertical wall" begin
            res = correct_collisions(
                [0, 0, 0], [0, 2, 0],
                mr.Walls(rotation=:x, position=[1])
            )
            compare(res, [
                [0, 0, 0], [0, 2, 0],
            ])
        end
        @testset "Hitting two vertical walls" begin
            res = correct_collisions(
                [0, 0, 0], [6, 0, 12],
                mr.Walls(rotation=:x, position=[-1, 1])
            )
            compare(res, [
                [0, 0, 0], [1, 0, 2], [-1 ,0, 6], [1, 0, 10], [0, 0, 12]
            ])
        end
        @testset "Hitting vertical and horizontal walls" begin
            res = correct_collisions(
                [-1, 0, 0], [2, 3, 3],
                [mr.Walls(rotation=:x, position=1.), mr.Walls(rotation=:y, position=1.)]
            )
            compare(res, [
                [-1, 0, 0], [0, 1, 1], [1, 0, 2], [0, -1, 3]
            ])
        end
        @testset "Hitting a corner" begin
            res = correct_collisions(
                [0, 0, 0], [3, 3, 3],
                [mr.Walls(rotation=:x, position=1.), mr.Walls(rotation=:y, position=1.)]
            )
            compare(res, [
                [0, 0, 0], [1, 1, 1], [1, 1, 1], [-1, -1, 3]
            ])
        end
        @testset "Hitting diagonal wall" begin
            res = correct_collisions(
                [1, 0, 0], [1, 3, 0],
                mr.Walls(rotation=[1, 1, 0], position=sqrt(2)),
            )
            compare(res, [
                [1, 0, 0], [1, 1, 0], [-1, 1, 0]
            ])
        end
        @testset "Test repeating wall" begin
            res = correct_collisions(
                [1.5, 0, 0], [3.5, 2, 0],
                mr.Walls(position=0., repeats=1.),
            )
            compare(res, [
                [1.5, 0, 0], [2, 0.5, 0], [1, 1.5, 0], [1.5, 2, 0]
            ])
        end
        @testset "Test shifted repeating wall" begin
            res = correct_collisions(
                [1, 0, 0], [5, 4, 0],
                mr.Walls(position=0.5, repeats=1.),
            )
            compare(res, [
                [1, 0, 0], 
                [1.5, 0.5, 0], [0.5, 1.5, 0], 
                [1.5, 2.5, 0], [0.5, 3.5, 0], 
                [1, 4., 0.],
            ])
        end
        @testset "Particle remain between reflecting walls" begin
            for (position, repeats) in [
                (0.5, 1.),
                (0.8, 1.6),
                (0., 11.3),
                (0., 10.7),
                (0.1, 1.6),
                (0.1, 11.3),
            ]
                @testset "$(repeats) um repeating wall placed at $position um" begin
                    Random.seed!(1234)
                    geometry = mr.Walls(position=position, repeats=repeats)
                    simulation = mr.Simulation([], geometry=geometry, diffusivity=3.)
                    snap = mr.Snapshot(10000)
                    snap2 = mr.evolve(snap, simulation, 3)
                    
                    get_pos(spin) = div(spin.position[1] - position, repeats, RoundDown)
                    @test all(get_pos.(snap) .== get_pos.(snap2))
                end
            end
        end
    end
    @testset "Sphere reflections" begin
        @testset "Remain within sphere" begin
            res = correct_collisions(
                [0, 0, 1.5], [6, 8, 6],
                mr.Spheres(radius=1., position=[0, 0, 2])
            )
            final = res[end]
            radius = norm(final .- [0, 0, 2])
            @assert radius <= 1.
        end
    end
    @testset "Cylinder reflections" begin
        @testset "Within cylinder along radial line" begin
            res = correct_collisions(
                [0, 0, 0], [6, 0, 6],
                [mr.Cylinders(radius=1., rotation=:z)]
            )
            compare(res, [
                [0, 0, 0], [1, 0, 1], [-1, 0, 3], [1, 0, 5], [0, 0, 6],
            ])
        end
        @testset "90 degree bounces within vertical cylinder" begin
            res = correct_collisions(
                [0, 1, 0], [10, 1, 0],
                mr.Cylinders(radius=sqrt(2), rotation=:z, position=[0, 0])
            )
            compare(res, [
                [0, 1, 0], [1, 1, 0], [1, -1, 0], [-1, -1, 0], [-1, 1, 0], [1, 1, 0], [1, 0, 0]
            ])
        end
        @testset "90 degree bounces from outside vertical cylinder" begin
            res = correct_collisions(
                [-2, 1, 0], [2, 1, 0],
                [mr.Cylinders(radius=sqrt(2), rotation=:z, position=[0, 0])]
            )
            compare(res, [
                [-2, 1, 0], [-1, 1, 0], [-1, 4, 0]
            ])
        end
        @testset "Remain within angled cylinder" begin
            orient = [1, 2, sqrt(3)]
            cylinder = mr.Cylinders(radius=2.3, rotation=orient)
            res = correct_collisions(
                [0, 0.5, 0.3], [-30, 50, 10],
                [cylinder]
            )
            final = res[end]
            radius = norm(final .- (orient ⋅ final) * orient / norm(orient) ^ 2)
            @test radius <= 2.3
            @test mr.isinside(cylinder, final) == 1
        end
        @testset "Remain within distant cylinder" begin
            Random.seed!(1234)
            geometry = mr.Cylinders(radius=[0.8, 0.9], repeats=[2., 2.])
            spin = mr.Spin(position=[200., 200., 0.])
            @test mr.isinside(geometry, spin) == 2
            inside = true
            empty_sequence = build_sequence() do 
                Sequence([10.]) 
            end
            seq_part = mr.SequenceParts.MultSequencePart(0.5, [mr.SequenceParts.EmptyPart()])
            for _ in 1:100
                mr.Evolve.draw_step!(spin, mr.Simulation(empty_sequence, diffusivity=3., geometry=geometry), seq_part, SizedVector{1, Float64}(3.))
                inside &= mr.isinside(geometry, spin) == 2
            end
            @test inside
        end
        @testset "Test variance of parallel displacement within cylinder" begin
            Random.seed!(1234)
            sim = mr.Simulation([], geometry=mr.Cylinders(radius=0.8), diffusivity=3.)
            snap = mr.evolve(zeros(10000, 3), sim, 10)
            zval = map(spin -> spin.position[3], snap)
            @test var(zval) ≈ 60. rtol=0.05
            @test maximum(spin -> norm(spin.position[1:2]), snap) ≈ 0.8 rtol=0.01
        end
    end
    @testset "Ray-grid intersections with undefined grid" begin
        function tcompare(t1, t2)
            @test length(t1) == length(t2)
            for (e1, e2) in zip(t1, t2)
                @test e1 ≈ e2
            end
        end
        res = collect(mr.Geometries.Internal.ray_grid_intersections([0.5, 0.5, 0.5], [0.5, 0.5, 3.5]))
        tcompare(res[1], ([0, 0, 0], 0., [0.5, 0.5, 0.5], 1/6, [0.5, 0.5, 1.]))
        tcompare(res[2], ([0, 0, 1], 1/6, [0.5, 0.5, 0.], 1/2, [0.5, 0.5, 1.]))
        tcompare(res[3], ([0, 0, 2], 1/2, [0.5, 0.5, 0.], 5/6, [0.5, 0.5, 1.]))
        tcompare(res[4], ([0, 0, 3], 5/6, [0.5, 0.5, 0.], 1., [0.5, 0.5, 0.5]))
    end
    @testset "Reflections on planes of cylinders" begin
        @testset "Bounce between four cylinders" begin
            cylinders = mr.Cylinders(
                radius=sqrt(2), repeats=[3, 4]
            )
            res = correct_collisions(
                [1, 2, 0], [1, 11, 9],
                [cylinders]
            )
            compare(res, [
                [1, 2, 0], [1, 3, 1], [2, 3, 2], [2, 1, 4],
                [1, 1, 5], [1, 3, 7], [2, 3, 8], [2, 2, 9],
            ])
        end
        @testset "Travel through many repeats between bounces" begin
            cylinders = mr.Cylinders(
                radius=1, repeats=[2, 4]
            )
            res = correct_collisions(
                [1, 2, 0], [13, 6, 4],
                [cylinders]
            )
            compare(res, [
                [1, 2, 0], [4, 3, 1], [10, 1, 3], [13, 2, 4]
            ])
        end
        @testset "Test lots of particles still don't cross compartments" begin
            for geometry in (
                mr.Cylinders(radius=[0.8, 0.9], position=[[0, 0], [0, 0]], repeats=[2, 2]),
                mr.Spheres(radius=[0.8, 0.9], position=[[0, 0, 0], [2, 0, 2]], repeats=[2, 2, 2]),
            )
                sequence = DWI(TE=80., bval=2., gradient=(type=:instant, ), diffusion_time=40.)

                snap = mr.Snapshot(300);

                simulation = mr.Simulation([sequence], diffusivity=3., geometry=geometry);

                before = mr.isinside(geometry, snap) .> 0
                after = mr.isinside(geometry, mr.evolve(snap, simulation, 200.)) .> 0
                switched = sum(xor.(before, after))
                @test switched == 0
            end
        end
        @testset "Test that we remain between two walls" begin
            Random.seed!(1234)
            walls = mr.Walls(position=[0, 1])
            snap = mr.Snapshot([mr.Spin(position=rand(3)) for _ in 1:3000])
            sequence = DWI(TE=80., bval=2., gradient=(type=:instant, ), diffusion_time=40.)
            simulation = mr.Simulation([sequence]; geometry=walls, diffusivity=3.);

            final = mr.evolve(snap, simulation, 200.)
            xfinal = [s.position[1] for s in final.spins]
            @test all(xfinal .>= 0.)
            @test all(xfinal .<= 1.)
        end
    end
    if false
        @testset "Test spiral collision detection" begin
            spiral = mr.spirals(1., 10., thickness=1., inner_cylinder=false)

            function test(origin, destination, distance, normal)
                c = mr.detect_collision(origin, destination, spiral, mr.empty_collision)
                @test c.distance ≈ distance
                @test (c.normal ./ norm(c.normal)) ≈ (SVector{3, Float64}(normal) ./ norm(normal))
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
            c = mr.detect_collision([0, 1, 0], [4, -1, 0], spiral, mr.empty_collision)
            @test 0 < c.distance < 0.25

            test([11, 1, 0], [7, -1, 0], 0.5, [1, 0, 0])
            c = mr.detect_collision([11, -1, 0], [7, 1, 0], spiral, mr.empty_collision)
            @test 0.25 < c.distance < 0.3
        end
        @testset "Spiral leakage detection" begin
            spiral = mr.spirals(0.8, 1., inner_cylinder=true, outer_cylinder=true)
            spin = mr.Spin(position=[0.9, 0., 0.])
            theta = mr.spiral_theta(spiral.obstructions[1], SVector{2}(spin.position[1:2]))
            for _ in 1:100000
                prev_theta = theta
                mr.Evolve.draw_step!(spin, 1., 0.01, (spiral, ), ())
                theta = mr.spiral_theta(spiral.obstructions[1], SVector{2}(spin.position[1:2]))
                if abs(theta - prev_theta) > π
                    error("Leaked through spiral!")
                end
                if mr.isinside(spiral, spin) != 1
                    error("Leaked through cylinders!")
                end
            end
        end
    end
end