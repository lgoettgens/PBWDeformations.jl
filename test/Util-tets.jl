@testset ExtendedTestSet "All Util.jl tests" begin
    @testset "groupBy" begin
        groupBy = PD.groupBy
        @test groupBy([1,1,2,1]) == [[1,1],[2],[1]]
    end
       
end
