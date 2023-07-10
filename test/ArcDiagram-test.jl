@testset "ArcDiagram.jl tests" begin
    @testset "numeric constructor" begin
        # wrong length
        @test_throws ArgumentError ArcDiagramUndirected(5, 5, [-2, -1, -4, -3, 1], [-5, 3, 2])
        @test_throws ArgumentError ArcDiagramUndirected(5, 5, [-2, -1, -4, -3], [1, -5, 3, 2, 5])
        @test_throws ArgumentError ArcDiagramUndirected(5, 5, [-2, -1, -4, -3], [1, -5, 3, 2, 5, 4])
        @test_throws ArgumentError ArcDiagramUndirected(5, 5, [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4, -1])

        # out of bounds
        @test_throws ArgumentError ArcDiagramUndirected(5, 5, [-6, -1, -4, -3, 2], [-5, 3, 2, 5, 4])
        @test_throws ArgumentError ArcDiagramUndirected(5, 5, [7, -1, -4, -3, 2], [-5, 3, 2, 5, 4])

        # self-loops
        @test_throws ArgumentError ArcDiagramUndirected(5, 5, [-1, -2, -4, -3, 1], [-5, 3, 2, 5, 4])
        @test_throws ArgumentError ArcDiagramUndirected(5, 5, [5, -1, -2, -3, -4], [-5, 1, 2, 3, 4])

        # correct
        @test ArcDiagramUndirected(5, 5, [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4]) !== nothing
    end

    @testset "string constructor" begin
        # not all symbols twice
        @test_throws ArgumentError ArcDiagramUndirected("ABC,ABD")
        @test_throws ArgumentError ArcDiagramUndirected("ABCD,ABDD")
        @test_throws ArgumentError ArcDiagramUndirected("ABDD,ABDD")

        # correct
        @test ArcDiagramUndirected("AACCE,EGGII") ==
              ArcDiagramUndirected("AACCE", "EGGII") ==
              ArcDiagramUndirected(5, 5, [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4])
        @test ArcDiagramUndirected("ABCD,ABCD") ==
              ArcDiagramUndirected("ABCD", "ABCD") ==
              ArcDiagramUndirected(4, 4, [1, 2, 3, 4], [-1, -2, -3, -4])
        @test ArcDiagramUndirected("ABBD,AD") ==
              ArcDiagramUndirected("ABBD", "AD") ==
              ArcDiagramUndirected(4, 2, [1, -3, -2, 2], [-1, -4])
    end

    @testset "is_crossing_free" begin
        @testset "is_crossing_free(part=:everywhere)" begin
            @test is_crossing_free(ArcDiagramUndirected("AA,CCEE")) == true
            @test is_crossing_free(ArcDiagramUndirected("AA,CDCD")) == false
            @test is_crossing_free(ArcDiagramUndirected("AA,CDDC")) == true
            @test is_crossing_free(ArcDiagramUndirected("AB,ABEE")) == true
            @test is_crossing_free(ArcDiagramUndirected("AB,ADBD")) == false
            @test is_crossing_free(ArcDiagramUndirected("AB,ADDB")) == true
            @test is_crossing_free(ArcDiagramUndirected("AB,BAEE")) == false
            @test is_crossing_free(ArcDiagramUndirected("AB,CABC")) == false
            @test is_crossing_free(ArcDiagramUndirected("AB,CACB")) == false
            @test is_crossing_free(ArcDiagramUndirected("AB,BDAD")) == false
            @test is_crossing_free(ArcDiagramUndirected("AB,CBAC")) == false
            @test is_crossing_free(ArcDiagramUndirected("AB,CCAB")) == true
            @test is_crossing_free(ArcDiagramUndirected("AB,BDDA")) == false
            @test is_crossing_free(ArcDiagramUndirected("AB,CBCA")) == false
            @test is_crossing_free(ArcDiagramUndirected("AB,CCBA")) == false

            @test is_crossing_free(ArcDiagramUndirected("AAC,CEE")) == true
            @test is_crossing_free(ArcDiagramUndirected("AAC,DCD")) == false
            @test is_crossing_free(ArcDiagramUndirected("AAC,DDC")) == true
            @test is_crossing_free(ArcDiagramUndirected("ABA,BEE")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABA,DBD")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABA,DDB")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABB,AEE")) == true
            @test is_crossing_free(ArcDiagramUndirected("ABC,ABC")) == true
            @test is_crossing_free(ArcDiagramUndirected("ABC,ACB")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABB,DAD")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABC,BAC")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABC,CAB")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABB,DDA")) == true
            @test is_crossing_free(ArcDiagramUndirected("ABC,BCA")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABC,CBA")) == false

            @test is_crossing_free(ArcDiagramUndirected("AACC,EE")) == true
            @test is_crossing_free(ArcDiagramUndirected("AACD,CD")) == true
            @test is_crossing_free(ArcDiagramUndirected("AACD,DC")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABAB,EE")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABAD,BD")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABAD,DB")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABBA,EE")) == true
            @test is_crossing_free(ArcDiagramUndirected("ABCA,BC")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABCA,CB")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABBD,AD")) == true
            @test is_crossing_free(ArcDiagramUndirected("ABCB,AC")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABCC,AB")) == true
            @test is_crossing_free(ArcDiagramUndirected("ABBD,DA")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABCB,CA")) == false
            @test is_crossing_free(ArcDiagramUndirected("ABCC,BA")) == false

        end

        @testset "is_crossing_free(part=:upper)" begin
            @test all(diag -> is_crossing_free(diag, part=:upper), all_arc_diagrams(0, 6))
            @test all(diag -> is_crossing_free(diag, part=:upper), all_arc_diagrams(1, 5))
            @test all(diag -> is_crossing_free(diag, part=:upper), all_arc_diagrams(2, 4))
            @test all(diag -> is_crossing_free(diag, part=:upper), all_arc_diagrams(3, 3))

            @test [diag for diag in all_arc_diagrams(4, 2) if !is_crossing_free(diag, part=:upper)] == [ArcDiagramUndirected("ABAB,EE")]

            @test is_crossing_free(ArcDiagramUndirected("AACCE,E"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("AACDC,D"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("AACDD,C"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("ABABE,E"), part=:upper) == false
            @test is_crossing_free(ArcDiagramUndirected("ABADB,D"), part=:upper) == false
            @test is_crossing_free(ArcDiagramUndirected("ABADD,B"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("ABBAE,E"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("ABCAB,C"), part=:upper) == false
            @test is_crossing_free(ArcDiagramUndirected("ABCAC,B"), part=:upper) == false
            @test is_crossing_free(ArcDiagramUndirected("ABBDA,D"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("ABCBA,C"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("ABCCA,B"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("ABBDD,A"), part=:upper) == true
            @test is_crossing_free(ArcDiagramUndirected("ABCBC,A"), part=:upper) == false
            @test is_crossing_free(ArcDiagramUndirected("ABCCB,A"), part=:upper) == true

            for diag in all_arc_diagrams(6, 0)
                @test is_crossing_free(diag, part=:upper) == is_crossing_free(diag)
            end
        end

        @testset "is_crossing_free(part=:lower)" begin
            @test all(diag -> is_crossing_free(diag, part=:lower), all_arc_diagrams(6, 0))
            @test all(diag -> is_crossing_free(diag, part=:lower), all_arc_diagrams(5, 1))
            @test all(diag -> is_crossing_free(diag, part=:lower), all_arc_diagrams(4, 2))
            @test all(diag -> is_crossing_free(diag, part=:lower), all_arc_diagrams(3, 3))

            @test [diag for diag in all_arc_diagrams(2, 4) if !is_crossing_free(diag, part=:lower)] == [ArcDiagramUndirected("AA,CDCD")]

            @test is_crossing_free(ArcDiagramUndirected("A,ACCEE"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,ACDCD"), part=:lower) == false
            @test is_crossing_free(ArcDiagramUndirected("A,ACDDC"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,BABEE"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,BADBD"), part=:lower) == false
            @test is_crossing_free(ArcDiagramUndirected("A,BADDB"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,BBAEE"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,BCABC"), part=:lower) == false
            @test is_crossing_free(ArcDiagramUndirected("A,BCACB"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,BBDAD"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,BCBAC"), part=:lower) == false
            @test is_crossing_free(ArcDiagramUndirected("A,BCCAB"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,BBDDA"), part=:lower) == true
            @test is_crossing_free(ArcDiagramUndirected("A,BCBCA"), part=:lower) == false
            @test is_crossing_free(ArcDiagramUndirected("A,BCCBA"), part=:lower) == true

            for diag in all_arc_diagrams(0, 6)
                @test is_crossing_free(diag, part=:lower) == is_crossing_free(diag)
            end
        end
    end

    @testset "all_arc_diagrams" begin
        for n in 0:5
            for i in 0:2n
                @test length(collect(all_arc_diagrams(i, 2n - i))) ==
                      length(all_arc_diagrams(i, 2n - i)) ==
                      div(factorial(2 * n), 2^n * factorial(n))
            end
            for i in 0:2n+1
                @test length(collect(all_arc_diagrams(i, 2n + 1 - i))) == length(all_arc_diagrams(i, 2n + 1 - i)) == 0
            end

            for i in 0:2n
                # two disjoint independent sets of size n each -> matchings on a bipartite graph
                indep_sets = [[-1:-1:-min(i, n); 1:n-min(i, n)], [-n-1:-1:-max(i, n); n-min(i, n)+1:2n-i]]
                @test length(collect(all_arc_diagrams(i, 2n - i; indep_sets))) ==
                      length(all_arc_diagrams(i, 2n - i; indep_sets)) ==
                      factorial(n)
            end
        end

        @test length(
                  collect(all_arc_diagrams(6, 6, indep_sets=[[-1, -2, -3], [-4, -5, -6], [1, 2], [3, 4], [5, 6]])),
              ) ==
              length(all_arc_diagrams(6, 6, indep_sets=[[-1, -2, -3], [-4, -5, -6], [1, 2], [3, 4], [5, 6]])) ==
              4440
    end

end
