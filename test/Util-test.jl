@testset "Util.jl tests" begin
    @testset "flatten" begin
        flatten = PD.flatten

        @test flatten([[]]) == []
        @test flatten(Vector{Any}[]) == []
        @test flatten([[1], [], [2, 3, 4], [5], []]) == 1:5
        @test flatten([[i] for i in 1:10]) == 1:10
        @test flatten([[i, i + 1] for i in 1:2:10]) == 1:10
        @test flatten([[i, i + 1, i + 2] for i in 1:3:10]) == 1:12
        @test flatten([collect(1:10)]) == 1:10
    end
end
