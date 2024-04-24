@testset "test_cli.jl" begin
import MCMRSimulator.CLI: run_main_test
import JSON
import CSV 
import DataFrames: DataFrame

function in_tmpdir(f)
    curdir = pwd()
    mktempdir() do fn
        cd(fn)
        try
            f()
        finally
            cd(curdir)
        end
    end
end

@testset "Test the creation commands do not crash" begin
    @testset "mcmr geometry create" begin
        @testset "create cylinders" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create cylinders 2 test.json --radius 1,2 --position 0,1:2,3 --R1_surface 1.5")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Cylinders"
                @test result["number"] == 2
                @test result["radius"] == [1., 2.]
                @test result["position"] == [[0., 1.], [2., 3.]]
                @test result["R1_surface"] == 1.5
                @test result["R1_inside"] == 0.
            end
        end
        @testset "create spheres" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create spheres 2 test.json --radius 1,2 --position 0,1,2:3,4,5 --R1_surface 1.5")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Spheres"
                @test result["number"] == 2
                @test result["radius"] == [1., 2.]
                @test result["position"] == [[0., 1., 2.], [3., 4., 5.]]
                @test result["R1_surface"] == 1.5
                @test result["R1_inside"] == 0.
            end
        end
        @testset "create walls" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create walls 2 test.json --position 0,1 --R1 1.5")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Walls"
                @test result["number"] == 2
                @test result["position"] == [0., 1.]
                @test result["R1"] == 1.5
                @test result["R2"] == 0.
            end
        end
        @testset "create annuli" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create annuli 2 test.json --inner 1,2 --outer 3 --myelin")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Annuli"
                @test result["number"] == 2
                @test result["position"] == [0., 0.]
                @test result["inner"] == [1, 2]
                @test result["outer"] == 3.
                @test result["myelin"] == true
            end
        end
        @testset "create bendy cylinder" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create bendy-cylinder 2 test.json --control_point 0,0,0:0,0.5,1 --radius 0.4,0.6 --closed 0,0,1 --repeats 1.3,2,2")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "BendyCylinder"
                @test result["number"] == 2
                @test result["control_point"] == [[0., 0., 0.], [0., 0.5, 1.]]
                @test result["radius"] == [0.4, 0.6]
                @test result["closed"] == [0, 0, 1]
                @test result["repeats"] == [1.3, 2., 2.]
            end
        end
        @testset "create based on file" begin
            in_tmpdir() do 
                open("positions.txt", "w") do f
                    write(f, "0 0
                    1 1
                    2 0.5")
                end
                open("radii.txt", "w") do f
                    write(f, "1.2 0.8 0.3")
                end
                open("R1.txt", "w") do f
                    write(f, "0.1")
                end
                _, err = run_main_test("geometry create cylinders 3 test.json --position positions.txt --radius radii.txt --R1_inside R1.txt")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Cylinders"
                @test result["number"] == 3
                @test result["position"] == [[0., 0.], [1., 1.], [2, 0.5]]
                @test result["radius"] == [1.2, 0.8, 0.3]
                @test result["R1_inside"] == 0.1
                @test result["R1_surface"] == 0.
            end
        end
    end
    @testset "mcmr geometry create-random" begin
        @testset "random cylinders" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random cylinders 0.5 test.json --repeats 20,20")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Cylinders"
                @test all(result["radius"] .≈ 1.)
                @test result["number"] > 3
            end
        end
        @testset "random spheres" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random spheres 0.5 test.json --repeats 20,20,20")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Spheres"
                @test all(result["radius"] .≈ 1.)
                @test result["number"] > 3
            end
        end
        @testset "random annuli" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random annuli 0.5 test.json --repeats 20,20 --g-ratio=0.8")
                @test length(err) == 0
                result = JSON.parse(open("test.json", "r"))
                @test result["type"] == "Annuli"
                @test all(result["outer"] .≈ 1.)
                @test all(result["inner"] .≈ 0.8)
                @test result["number"] > 3
            end
        end
    end
end
@testset "Test full simulations" begin
    @testset "Setting R2" begin
        in_tmpdir() do
            @testset "Single global R2" begin
                _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2,2.2,2.2")
                @test length(err) == 0
                gradient_echo = GradientEcho(TE=30.)
                write_sequence("ge.seq", gradient_echo)
                @test length(err) == 0
                _, err = run_main_test("run spheres.json ge.seq --R2 0.1 -N 100 -o global.csv")
                @test length(err) == 0
                result = DataFrame(CSV.File("global.csv"))
                @test size(result, 1) == 1
                @test result[1, :nspins] == 100
                @test result[1, :transverse] ≈ 100 * exp(-3.)
            end
            @testset "Change R2 within sphere" begin
                _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2,2.2,2.2 --R2_inside=0.2")
                @test length(err) == 0
                _, err = run_main_test("run spheres.json ge.seq --R2 0.1 -N 100 -o varies.csv --subset inside --subset outside")
                @test length(err) == 0
                result = DataFrame(CSV.File("varies.csv"))
                @test size(result, 1) == 3
                @test result[1, :nspins] == 100
                @test result[2, :nspins] < 100
                @test result[3, :nspins] < 100
                @test result[2, :transverse] / result[2, :nspins] ≈ exp(-9.)
                @test result[3, :transverse] / result[3, :nspins] ≈ exp(-3.)
                @test all(result[!, :sequence] .== "ge.seq")
                @test all(result[!, :sequence_index] .== 1)
            end
        end
    end
end
@testset "Setting R1" begin
    in_tmpdir() do
        _, err = run_main_test("geometry create spheres 1 spheres.json --radius 1 --repeats 2.2,2.2,2.2 --R1_inside=0.02")
        @test length(err) == 0
        gradient_echo = GradientEcho(TE=30.)
        write_sequence("ge.seq", gradient_echo)
        _, err = run_main_test("run spheres.json ge.seq --R1 0.01 -N 100 -o R1.csv --subset inside --subset outside")
        @test length(err) == 0
        result = DataFrame(CSV.File("R1.csv"))
        @test size(result, 1) == 3
        @test result[1, :nspins] == 100
        @test result[2, :nspins] < 100
        @test result[3, :nspins] < 100
        @test result[2, :longitudinal] / result[2, :nspins] ≈ 1. - exp(-0.9)
        @test result[3, :longitudinal] / result[3, :nspins] ≈ 1. - exp(-0.3)
        @test all(result[!, :sequence] .== "ge.seq")
        @test all(result[!, :sequence_index] .== 1)
    end
end
if false
@testset "Using bvecs" begin
    in_tmpdir() do
        _, err = run_main_test("geometry create walls 1 walls.json --repeats 1")
        @test length(err) == 0
        write_sequence("dwi.seq", DWI(TE=80., bval=0.5))
        @test length(err) == 0
        open("bvecs", "w") do f
            write(f, "1 0
            0 1
            0 0")
        end
        _, err = run_main_test("run walls.json dwi.seq -o with_diff.csv --diffusivity 3. --bvecs=bvecs")
        @test length(err) == 0
        with_diff = DataFrame(CSV.File("with_diff.csv"))
        @test size(with_diff, 1) == 2
        @test with_diff[!, :sequence_index] == [1, 1]
        @test with_diff[!, :sequence] == ["dwi.json", "dwi.json"]
        @test with_diff[!, :bvec] == [1, 2]
        @test with_diff[1, :transverse] > exp(-1) * with_diff[1, :nspins]
        @test with_diff[2, :transverse] ≈ with_diff[2, :nspins] * exp(-1.5) rtol=0.1

        _, err = run_main_test("run walls.json dwi.json -o no_diff.csv --diffusivity 0. --bvecs=bvecs")
        no_diff = DataFrame(CSV.File("no_diff.csv"))
        @test size(no_diff, 1) == 2
        @test with_diff[!, :sequence_index] == [1, 1]
        @test with_diff[!, :sequence] == ["dwi.json", "dwi.json"]
        @test no_diff[!, :bvec] == [1, 2]
        @test no_diff[1, :transverse] ≈ no_diff[1, :nspins]
        @test no_diff[2, :transverse] ≈ no_diff[2, :nspins]
    end
end
end

end
