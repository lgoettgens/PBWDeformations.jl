@testset ExtendedTestSet "All PBWDeformations.Combinatorics tests" begin
    @testset "multichoose" begin
        for k in 1:10
            @test PD.multichoose(0, k) == 0
            @test PD.multichoose(1, k) == 1
            @test PD.multichoose(2, k) == k+1
        end
        for n in 1:10
            @test PD.multichoose(n, 0) == 1
            @test PD.multichoose(n, 1) == n
            @test PD.multichoose(n, 2) == div(n*(n+1), 2)
        end

        @test PD.multichoose(3, 3) == 10
    end

    @testset "multicombinations" begin
        for _ in 1:numRandomTests
            n = rand(0:20)
            k = rand(-1:10)

            @test PD.multichoose(n, k) ==
                length(PD.multicombinations(1:n, k)) ==
                length(PD.multicombinations(collect(1:n), k)) ==
                length(collect(PD.multicombinations(1:n, k))) ==
                length(collect(PD.multicombinations(collect(1:n), k)))
        end

        @test map(join, PD.multicombinations(["a", "b", "c"], 0)) == [""]
        @test map(join, PD.multicombinations(["a", "b", "c"], 1)) == ["a", "b", "c"]
        @test map(join, PD.multicombinations(["a", "b", "c"], 2)) == ["aa", "ab", "ac", "bb", "bc", "cc"]
        @test map(join, PD.multicombinations(["a", "b", "c"], 3)) == ["aaa", "aab", "aac", "abb", "abc", "acc", "bbb", "bbc", "bcc", "ccc"]
    end
    
end
