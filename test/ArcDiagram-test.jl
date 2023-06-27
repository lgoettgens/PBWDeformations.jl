@testset "ArcDiagram.jl tests" begin
    @testset "numeric constructor" begin
        # wrong length
        @test_throws ArgumentError ArcDiagram(5, 5, [-2, -1, -4, -3, 1], [-5, 3, 2])
        @test_throws ArgumentError ArcDiagram(5, 5, [-2, -1, -4, -3], [1, -5, 3, 2, 5])
        @test_throws ArgumentError ArcDiagram(5, 5, [-2, -1, -4, -3], [1, -5, 3, 2, 5, 4])
        @test_throws ArgumentError ArcDiagram(5, 5, [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4, -1])

        # out of bounds
        @test_throws ArgumentError ArcDiagram(5, 5, [-6, -1, -4, -3, 2], [-5, 3, 2, 5, 4])
        @test_throws ArgumentError ArcDiagram(5, 5, [7, -1, -4, -3, 2], [-5, 3, 2, 5, 4])

        # self-loops
        @test_throws ArgumentError ArcDiagram(5, 5, [-1, -2, -4, -3, 1], [-5, 3, 2, 5, 4])
        @test_throws ArgumentError ArcDiagram(5, 5, [5, -1, -2, -3, -4], [-5, 1, 2, 3, 4])

        # correct
        @test ArcDiagram(5, 5, [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4]) !== nothing
    end

    @testset "string constructor" begin
        # not all symbols twice
        @test_throws ArgumentError ArcDiagram("ABC,ABD")
        @test_throws ArgumentError ArcDiagram("ABCD,ABDD")
        @test_throws ArgumentError ArcDiagram("ABDD,ABDD")

        # correct
        @test ArcDiagram("AACCE,EGGII") ==
              ArcDiagram("AACCE", "EGGII") ==
              ArcDiagram(5, 5, [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4])
        @test ArcDiagram("ABCD,ABCD") == ArcDiagram("ABCD", "ABCD") == ArcDiagram(4, 4, [1, 2, 3, 4], [-1, -2, -3, -4])
        @test ArcDiagram("ABBD,AD") == ArcDiagram("ABBD", "AD") == ArcDiagram(4, 2, [1, -3, -2, 2], [-1, -4])
    end

    @testset "is_crossing_free" begin
        @testset "is_crossing_free(part=:everywhere)" begin
            @test is_crossing_free(ArcDiagram("AA,CCEE")) == true
            @test is_crossing_free(ArcDiagram("AA,CDCD")) == false
            @test is_crossing_free(ArcDiagram("AA,CDDC")) == true
            @test is_crossing_free(ArcDiagram("AB,ABEE")) == true
            @test is_crossing_free(ArcDiagram("AB,ADBD")) == false
            @test is_crossing_free(ArcDiagram("AB,ADDB")) == true
            @test is_crossing_free(ArcDiagram("AB,BAEE")) == false
            @test is_crossing_free(ArcDiagram("AB,CABC")) == false
            @test is_crossing_free(ArcDiagram("AB,CACB")) == false
            @test is_crossing_free(ArcDiagram("AB,BDAD")) == false
            @test is_crossing_free(ArcDiagram("AB,CBAC")) == false
            @test is_crossing_free(ArcDiagram("AB,CCAB")) == true
            @test is_crossing_free(ArcDiagram("AB,BDDA")) == false
            @test is_crossing_free(ArcDiagram("AB,CBCA")) == false
            @test is_crossing_free(ArcDiagram("AB,CCBA")) == false

            @test is_crossing_free(ArcDiagram("AAC,CEE")) == true
            @test is_crossing_free(ArcDiagram("AAC,DCD")) == false
            @test is_crossing_free(ArcDiagram("AAC,DDC")) == true
            @test is_crossing_free(ArcDiagram("ABA,BEE")) == false
            @test is_crossing_free(ArcDiagram("ABA,DBD")) == false
            @test is_crossing_free(ArcDiagram("ABA,DDB")) == false
            @test is_crossing_free(ArcDiagram("ABB,AEE")) == true
            @test is_crossing_free(ArcDiagram("ABC,ABC")) == true
            @test is_crossing_free(ArcDiagram("ABC,ACB")) == false
            @test is_crossing_free(ArcDiagram("ABB,DAD")) == false
            @test is_crossing_free(ArcDiagram("ABC,BAC")) == false
            @test is_crossing_free(ArcDiagram("ABC,CAB")) == false
            @test is_crossing_free(ArcDiagram("ABB,DDA")) == true
            @test is_crossing_free(ArcDiagram("ABC,BCA")) == false
            @test is_crossing_free(ArcDiagram("ABC,CBA")) == false

            @test is_crossing_free(ArcDiagram("AACC,EE")) == true
            @test is_crossing_free(ArcDiagram("AACD,CD")) == true
            @test is_crossing_free(ArcDiagram("AACD,DC")) == false
            @test is_crossing_free(ArcDiagram("ABAB,EE")) == false
            @test is_crossing_free(ArcDiagram("ABAD,BD")) == false
            @test is_crossing_free(ArcDiagram("ABAD,DB")) == false
            @test is_crossing_free(ArcDiagram("ABBA,EE")) == true
            @test is_crossing_free(ArcDiagram("ABCA,BC")) == false
            @test is_crossing_free(ArcDiagram("ABCA,CB")) == false
            @test is_crossing_free(ArcDiagram("ABBD,AD")) == true
            @test is_crossing_free(ArcDiagram("ABCB,AC")) == false
            @test is_crossing_free(ArcDiagram("ABCC,AB")) == true
            @test is_crossing_free(ArcDiagram("ABBD,DA")) == false
            @test is_crossing_free(ArcDiagram("ABCB,CA")) == false
            @test is_crossing_free(ArcDiagram("ABCC,BA")) == false

        end

        @testset "is_crossing_free(part=:upper)" begin
            @test all(diag -> is_crossing_free(diag, part=:upper), all_arc_diagrams(0, 6))
            @test all(diag -> is_crossing_free(diag, part=:upper), all_arc_diagrams(1, 5))
            @test all(diag -> is_crossing_free(diag, part=:upper), all_arc_diagrams(2, 4))
            @test all(diag -> is_crossing_free(diag, part=:upper), all_arc_diagrams(3, 3))

            @test [diag for diag in all_arc_diagrams(4, 2) if !is_crossing_free(diag, part=:upper)] == [ArcDiagram("ABAB,EE")]

            @test is_crossing_free(ArcDiagram("AACCE,E"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("AACDC,D"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("AACDD,C"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("ABABE,E"), part=:upper) == false
            @test is_crossing_free(ArcDiagram("ABADB,D"), part=:upper) == false
            @test is_crossing_free(ArcDiagram("ABADD,B"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("ABBAE,E"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("ABCAB,C"), part=:upper) == false
            @test is_crossing_free(ArcDiagram("ABCAC,B"), part=:upper) == false
            @test is_crossing_free(ArcDiagram("ABBDA,D"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("ABCBA,C"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("ABCCA,B"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("ABBDD,A"), part=:upper) == true
            @test is_crossing_free(ArcDiagram("ABCBC,A"), part=:upper) == false
            @test is_crossing_free(ArcDiagram("ABCCB,A"), part=:upper) == true

            for diag in all_arc_diagrams(6, 0)
                @test is_crossing_free(diag, part=:upper) == is_crossing_free(diag)
            end
        end

        @testset "is_crossing_free(part=:lower)" begin
            @test all(diag -> is_crossing_free(diag, part=:lower), all_arc_diagrams(6, 0))
            @test all(diag -> is_crossing_free(diag, part=:lower), all_arc_diagrams(5, 1))
            @test all(diag -> is_crossing_free(diag, part=:lower), all_arc_diagrams(4, 2))
            @test all(diag -> is_crossing_free(diag, part=:lower), all_arc_diagrams(3, 3))

            @test [diag for diag in all_arc_diagrams(2, 4) if !is_crossing_free(diag, part=:lower)] == [ArcDiagram("AA,CDCD")]

            @test is_crossing_free(ArcDiagram("A,ACCEE"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,ACDCD"), part=:lower) == false
            @test is_crossing_free(ArcDiagram("A,ACDDC"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,BABEE"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,BADBD"), part=:lower) == false
            @test is_crossing_free(ArcDiagram("A,BADDB"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,BBAEE"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,BCABC"), part=:lower) == false
            @test is_crossing_free(ArcDiagram("A,BCACB"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,BBDAD"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,BCBAC"), part=:lower) == false
            @test is_crossing_free(ArcDiagram("A,BCCAB"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,BBDDA"), part=:lower) == true
            @test is_crossing_free(ArcDiagram("A,BCBCA"), part=:lower) == false
            @test is_crossing_free(ArcDiagram("A,BCCBA"), part=:lower) == true

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
