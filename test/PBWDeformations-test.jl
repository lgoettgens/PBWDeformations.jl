@testset ExtendedTestSet "All PBWDeformations.jl tests" begin
    @testset "isvalidlie_for_gap" begin

        function testit(dynkin, pred, until=10)
            for i in 0:until
                @test PD.isvalidlie_for_gap(dynkin, i) == pred(i)
            end
        end

        testit('A', >=(1))
        testit('B', >=(2))
        testit('C', >=(2))
        testit('D', >=(4))
        testit('E', in([6,7,8]))
        testit('F', ==(4))
        testit('G', ==(2))
    end
       
end
