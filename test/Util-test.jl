@testset ExtendedTestSet "All Util.jl tests" begin
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

    @testset "groupBy" begin
        groupBy = PD.groupBy

        @test groupBy([]) == []
        @test groupBy([1, 1, 2, 1]) == [[1, 1], [2], [1]]
        @test groupBy([1, 1, 2, 2, 2, 2, 3, 1, 4, 4]) == [[1, 1], [2, 2, 2, 2], [3], [1], [4, 4]]
        @test groupBy([1 for _ in 1:10]) == [[1 for _ in 1:10]]
        @test groupBy([i for i in 1:10]) == [[i] for i in 1:10]

        @test groupBy([i for i in -5:5]; eq=((x, y) -> sign(x) == sign(y))) == [[i for i in -5:-1], [0], [i for i in 1:5]]
    end

    @testset "is_valid_dynkin" begin
        function testit(dynkin, pred, until=10)
            for i in 0:until
                @test PD.is_valid_dynkin(dynkin, i) == pred(i)
            end
        end

        testit('A', >=(1))
        testit('B', >=(2))
        testit('C', >=(2))
        testit('D', >=(4))
        testit('E', in([6, 7, 8]))
        testit('F', ==(4))
        testit('G', ==(2))
        testit('X', _ -> false)
    end
end
