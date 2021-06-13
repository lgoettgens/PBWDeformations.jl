lie = PD.lie
mod = PD.mod

@testset ExtendedTestSet "All PBWDeformations.jl tests" begin
    @testset "sanitizeLieInput" begin
        @test_throws AssertionError PD.sanitizeLieInput('A', 0)
        @test_throws AssertionError PD.sanitizeLieInput('B', 0)
        @test_throws AssertionError PD.sanitizeLieInput('C', 0)
        @test_throws AssertionError PD.sanitizeLieInput('D', 0)

        @test_throws AssertionError PD.sanitizeLieInput('B', 1)
        @test_throws AssertionError PD.sanitizeLieInput('C', 1)
        @test_throws AssertionError PD.sanitizeLieInput('D', 1)

        @test_throws AssertionError PD.sanitizeLieInput('D', 2)
        @test_throws AssertionError PD.sanitizeLieInput('D', 3)

        @test_nowarn PD.sanitizeLieInput('A',1)
        @test_nowarn PD.sanitizeLieInput('B',2)
        @test_nowarn PD.sanitizeLieInput('C',2)
        @test_nowarn PD.sanitizeLieInput('D',4)
    end
       
end
