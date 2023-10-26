@testset "test_cli.jl" begin
import MCMRSimulator.CLI: run_main_test
import JSON: parse
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
                _, err = run_main_test("geometry create cylinders 2 test.json --radius 1 2 --position 0 1 2 3 --R1_surface 1.5")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
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
                _, err = run_main_test("geometry create spheres 2 test.json --radius 1 2 --position 0 1 2 3 4 5 --R1_surface 1.5")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
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
                _, err = run_main_test("geometry create walls 2 test.json --position 0 1 --R1 1.5")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Walls"
                @test result["number"] == 2
                @test result["position"] == [0., 1.]
                @test result["R1"] == 1.5
                @test result["R2"] == 0.
            end
        end
        @testset "create annuli" begin
            in_tmpdir() do 
                _, err = run_main_test("geometry create annuli 2 test.json --inner 1 2 --outer 3 --myelin")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Annuli"
                @test result["number"] == 2
                @test result["position"] == [0., 0.]
                @test result["inner"] == [1, 2]
                @test result["outer"] == 3.
                @test result["myelin"] == true
            end
        end
    end
    @testset "mcmr geometry create-random" begin
        @testset "random cylinders" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random cylinders 0.5 test.json --repeats 20 20")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Cylinders"
                @test all(result["radius"] .≈ 1.)
                @test result["number"] > 3
            end
        end
        @testset "random spheres" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random spheres 0.5 test.json --repeats 20 20 20")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Spheres"
                @test all(result["radius"] .≈ 1.)
                @test result["number"] > 3
            end
        end
        @testset "random annuli" begin
            in_tmpdir() do
                _, err = run_main_test("geometry create-random annuli 0.5 test.json --repeats 20 20 --g-ratio=0.8")
                @test length(err) == 0
                result = parse(open("test.json", "r"))
                @test result["type"] == "Annuli"
                @test all(result["outer"] .≈ 1.)
                @test all(result["inner"] .≈ 0.8)
                @test result["number"] > 3
            end
        end
    end
end

end
