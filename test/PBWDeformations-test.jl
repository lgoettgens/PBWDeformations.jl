@testset ExtendedTestSet "All PBWDeformations.jl tests" begin
    @testset "sanitize_lie_input" begin
        @test_throws AssertionError PD.sanitize_lie_input('A', 0)
        @test_throws AssertionError PD.sanitize_lie_input('B', 0)
        @test_throws AssertionError PD.sanitize_lie_input('C', 0)
        @test_throws AssertionError PD.sanitize_lie_input('D', 0)

        @test_throws AssertionError PD.sanitize_lie_input('B', 1)
        @test_throws AssertionError PD.sanitize_lie_input('C', 1)
        @test_throws AssertionError PD.sanitize_lie_input('D', 1)

        @test_throws AssertionError PD.sanitize_lie_input('D', 2)
        @test_throws AssertionError PD.sanitize_lie_input('D', 3)

        @test_nowarn PD.sanitize_lie_input('A',1)
        @test_nowarn PD.sanitize_lie_input('B',2)
        @test_nowarn PD.sanitize_lie_input('C',2)
        @test_nowarn PD.sanitize_lie_input('D',4)
    end
       
end
