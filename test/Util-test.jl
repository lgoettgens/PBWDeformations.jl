std_basis = PD.std_basis
ur_triag_entries = PD.ur_triag_entries
ur_proper_triag_entries = PD.ur_proper_triag_entries

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


    @testset "std_basis" begin
        @test std_basis(1, 3) == [1, 0, 0]
        @test std_basis(2, 3) == [0, 1, 0]
        @test std_basis(3, 3) == [0, 0, 1]

        @test std_basis(1, 4) == [1, 0, 0, 0]
        @test std_basis(2, 4) == [0, 1, 0, 0]
        @test std_basis(3, 4) == [0, 0, 1, 0]
        @test std_basis(4, 4) == [0, 0, 0, 1]

        for n in 1:10
            for i in 1:n
                @test sum(std_basis(i, n)) == 1
                for j in 1:n
                    @test std_basis(i, n)[j] == (i == j ? 1 : 0)
                end
            end
        end

        @testset "std_basis, T = $T" for T in [Int, Float64, UInt16, Bool, Rational{Int}, Complex{Int}]
            @test std_basis(T, 1, 3) == [one(T), zero(T), zero(T)]
            @test std_basis(T, 2, 3) == [zero(T), one(T), zero(T)]
            @test std_basis(T, 3, 3) == [zero(T), zero(T), one(T)]

            @test std_basis(T, 1, 4) == [one(T), zero(T), zero(T), zero(T)]
            @test std_basis(T, 2, 4) == [zero(T), one(T), zero(T), zero(T)]
            @test std_basis(T, 3, 4) == [zero(T), zero(T), one(T), zero(T)]
            @test std_basis(T, 4, 4) == [zero(T), zero(T), zero(T), one(T)]

            for n in 1:6
                for i in 1:n
                    @test sum(std_basis(T, i, n)) == one(T)
                    for j in 1:n
                        @test std_basis(T, i, n)[j] == (i == j ? one(T) : zero(T))
                    end
                end
            end
        end
    end


    @testset "ur_triag_entries" begin
        @test ur_triag_entries([1 2; 3 4]) == [1, 2, 4]
        @test ur_triag_entries([1 2 3; 4 5 6; 7 8 9]) == [1, 2, 3, 5, 6, 9]

        for n in 1:10
            M = zeros(n, n)
            @test length(ur_triag_entries(M)) == div(n * (n + 1), 2)
        end
    end


    @testset "ur_proper_triag_entries" begin
        @test ur_proper_triag_entries([1 2; 3 4]) == [2]
        @test ur_proper_triag_entries([1 2 3; 4 5 6; 7 8 9]) == [2, 3, 6]

        for n in 1:10
            M = zeros(n, n)
            @test length(ur_proper_triag_entries(M)) == div(n * (n - 1), 2)
        end
    end

end
