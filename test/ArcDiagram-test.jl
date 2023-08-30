@testset "ArcDiagram.jl tests" begin
    @testset "Undirected" begin
        @testset "numeric constructor" begin
            # out of bounds
            @test_throws ArgumentError arc_diagram(Undirected, [-6, -1, -4, -3, 2], [-5, 3, 2, 5, 4])
            @test_throws ArgumentError arc_diagram(Undirected, [7, -1, -4, -3, 2], [-5, 3, 2, 5, 4])

            # self-loops
            @test_throws ArgumentError arc_diagram(Undirected, [-1, -2, -4, -3, 1], [-5, 3, 2, 5, 4])
            @test_throws ArgumentError arc_diagram(Undirected, [5, -1, -2, -3, -4], [-5, 1, 2, 3, 4])

            # correct
            @test arc_diagram(Undirected, [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4]) !== nothing
        end

        @testset "string constructor" begin
            # not all symbols twice
            @test_throws ArgumentError arc_diagram(Undirected, "ABC,ABD")
            @test_throws ArgumentError arc_diagram(Undirected, "ABCD,ABDD")
            @test_throws ArgumentError arc_diagram(Undirected, "ABDD,ABDD")

            # correct
            @test arc_diagram(Undirected, "AACCE,EGGII") ==
                  arc_diagram(Undirected, "AACCE", "EGGII") ==
                  arc_diagram(Undirected, [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4])
            @test arc_diagram(Undirected, "ABCD,ABCD") ==
                  arc_diagram(Undirected, "ABCD", "ABCD") ==
                  arc_diagram(Undirected, [1, 2, 3, 4], [-1, -2, -3, -4])
            @test arc_diagram(Undirected, "ABBD,AD") ==
                  arc_diagram(Undirected, "ABBD", "AD") ==
                  arc_diagram(Undirected, [1, -3, -2, 2], [-1, -4])
        end

        @testset "is_crossing_free" begin
            @testset "is_crossing_free(part=:everywhere)" begin
                @test is_crossing_free(arc_diagram(Undirected, "AA,CCEE")) == true
                @test is_crossing_free(arc_diagram(Undirected, "AA,CDCD")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AA,CDDC")) == true
                @test is_crossing_free(arc_diagram(Undirected, "AB,ABEE")) == true
                @test is_crossing_free(arc_diagram(Undirected, "AB,ADBD")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AB,ADDB")) == true
                @test is_crossing_free(arc_diagram(Undirected, "AB,BAEE")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AB,CABC")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AB,CACB")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AB,BDAD")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AB,CBAC")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AB,CCAB")) == true
                @test is_crossing_free(arc_diagram(Undirected, "AB,BDDA")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AB,CBCA")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AB,CCBA")) == false

                @test is_crossing_free(arc_diagram(Undirected, "AAC,CEE")) == true
                @test is_crossing_free(arc_diagram(Undirected, "AAC,DCD")) == false
                @test is_crossing_free(arc_diagram(Undirected, "AAC,DDC")) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABA,BEE")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABA,DBD")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABA,DDB")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABB,AEE")) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABC,ABC")) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABC,ACB")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABB,DAD")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABC,BAC")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABC,CAB")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABB,DDA")) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABC,BCA")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABC,CBA")) == false

                @test is_crossing_free(arc_diagram(Undirected, "AACC,EE")) == true
                @test is_crossing_free(arc_diagram(Undirected, "AACD,CD")) == true
                @test is_crossing_free(arc_diagram(Undirected, "AACD,DC")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABAB,EE")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABAD,BD")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABAD,DB")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABBA,EE")) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABCA,BC")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABCA,CB")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABBD,AD")) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABCB,AC")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABCC,AB")) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABBD,DA")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABCB,CA")) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABCC,BA")) == false

            end

            @testset "is_crossing_free(part=:upper)" begin
                @test all(diag -> is_crossing_free(diag; part=:upper), all_arc_diagrams(Undirected, 0, 6))
                @test all(diag -> is_crossing_free(diag; part=:upper), all_arc_diagrams(Undirected, 1, 5))
                @test all(diag -> is_crossing_free(diag; part=:upper), all_arc_diagrams(Undirected, 2, 4))
                @test all(diag -> is_crossing_free(diag; part=:upper), all_arc_diagrams(Undirected, 3, 3))

                @test [diag for diag in all_arc_diagrams(Undirected, 4, 2) if !is_crossing_free(diag; part=:upper)] == [arc_diagram(Undirected, "ABAB,EE")]

                @test is_crossing_free(arc_diagram(Undirected, "AACCE,E"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "AACDC,D"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "AACDD,C"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABABE,E"); part=:upper) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABADB,D"); part=:upper) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABADD,B"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABBAE,E"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABCAB,C"); part=:upper) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABCAC,B"); part=:upper) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABBDA,D"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABCBA,C"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABCCA,B"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABBDD,A"); part=:upper) == true
                @test is_crossing_free(arc_diagram(Undirected, "ABCBC,A"); part=:upper) == false
                @test is_crossing_free(arc_diagram(Undirected, "ABCCB,A"); part=:upper) == true

                for diag in all_arc_diagrams(Undirected, 6, 0)
                    @test is_crossing_free(diag; part=:upper) == is_crossing_free(diag)
                end
            end

            @testset "is_crossing_free(part=:lower)" begin
                @test all(diag -> is_crossing_free(diag; part=:lower), all_arc_diagrams(Undirected, 6, 0))
                @test all(diag -> is_crossing_free(diag; part=:lower), all_arc_diagrams(Undirected, 5, 1))
                @test all(diag -> is_crossing_free(diag; part=:lower), all_arc_diagrams(Undirected, 4, 2))
                @test all(diag -> is_crossing_free(diag; part=:lower), all_arc_diagrams(Undirected, 3, 3))

                @test [diag for diag in all_arc_diagrams(Undirected, 2, 4) if !is_crossing_free(diag; part=:lower)] == [arc_diagram(Undirected, "AA,CDCD")]

                @test is_crossing_free(arc_diagram(Undirected, "A,ACCEE"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,ACDCD"); part=:lower) == false
                @test is_crossing_free(arc_diagram(Undirected, "A,ACDDC"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,BABEE"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,BADBD"); part=:lower) == false
                @test is_crossing_free(arc_diagram(Undirected, "A,BADDB"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,BBAEE"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,BCABC"); part=:lower) == false
                @test is_crossing_free(arc_diagram(Undirected, "A,BCACB"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,BBDAD"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,BCBAC"); part=:lower) == false
                @test is_crossing_free(arc_diagram(Undirected, "A,BCCAB"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,BBDDA"); part=:lower) == true
                @test is_crossing_free(arc_diagram(Undirected, "A,BCBCA"); part=:lower) == false
                @test is_crossing_free(arc_diagram(Undirected, "A,BCCBA"); part=:lower) == true

                for diag in all_arc_diagrams(Undirected, 0, 6)
                    @test is_crossing_free(diag; part=:lower) == is_crossing_free(diag)
                end
            end
        end

        @testset "all_arc_diagrams(Undirected, ...)" begin
            for n in 0:5
                for i in 0:2n
                    @test length(collect(all_arc_diagrams(Undirected, i, 2n - i))) ==
                          length(all_arc_diagrams(Undirected, i, 2n - i)) ==
                          div(factorial(2 * n), 2^n * factorial(n))
                end
                for i in 0:2n+1
                    @test length(collect(all_arc_diagrams(Undirected, i, 2n + 1 - i))) ==
                          length(all_arc_diagrams(Undirected, i, 2n + 1 - i)) ==
                          0
                end

                for i in 0:2n
                    # two disjoint independent sets of size n each -> matchings on a bipartite graph
                    indep_sets = [[-1:-1:-min(i, n); 1:n-min(i, n)], [-n-1:-1:-max(i, n); n-min(i, n)+1:2n-i]]
                    @test length(collect(all_arc_diagrams(Undirected, i, 2n - i; indep_sets))) ==
                          length(all_arc_diagrams(Undirected, i, 2n - i; indep_sets)) ==
                          factorial(n)
                end
            end

            @test length(
                      all_arc_diagrams(
                          Undirected,
                          6,
                          6,
                          indep_sets=[[-1, -2, -3], [-4, -5, -6], [1, 2], [3, 4], [5, 6]],
                      ),
                  ) ==
                  length(
                      collect(
                          all_arc_diagrams(
                              Undirected,
                              6,
                              6,
                              indep_sets=[[-1, -2, -3], [-4, -5, -6], [1, 2], [3, 4], [5, 6]],
                          ),
                      ),
                  ) ==
                  4440
        end

    end

    @testset "Directed" begin
        @testset "numeric constructor" begin
            # out of bounds
            @test_throws ArgumentError arc_diagram(
                Directed,
                [1, 0, 1, 0, 1],
                [0, 1, 0, 1, 0],
                [-6, -1, -4, -3, 2],
                [-5, 3, 2, 5, 4],
            )
            @test_throws ArgumentError arc_diagram(
                Directed,
                [1, 0, 1, 0, 1],
                [0, 1, 0, 1, 0],
                [7, -1, -4, -3, 2],
                [-5, 3, 2, 5, 4],
            )

            # self-loops
            @test_throws ArgumentError arc_diagram(
                Directed,
                [1, 0, 1, 0, 1],
                [0, 1, 0, 1, 0],
                [-1, -2, -4, -3, 1],
                [-5, 3, 2, 5, 4],
            )
            @test_throws ArgumentError arc_diagram(
                Directed,
                [1, 0, 1, 0, 1],
                [0, 1, 0, 1, 0],
                [5, -1, -2, -3, -4],
                [-5, 1, 2, 3, 4],
            )

            # parity mismatch
            @test_throws ArgumentError arc_diagram(
                Directed,
                [1, 0, 1, 0, 1],
                [0, 1, 0, 1, 0],
                [-2, -1, -4, -3, 1],
                [-5, 3, 2, 5, 4],
            )

            @test_throws ArgumentError arc_diagram(
                Directed,
                [1, 1, 1, 1, 1],
                [0, 0, 0, 1, 0],
                [-2, -1, -4, -3, 1],
                [-5, 3, 2, 5, 4],
            )

            # correct
            @test arc_diagram(Directed, [1, 0, 0, 1, 1], [1, 1, 0, 0, 1], [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4]) !==
                  nothing
        end

        @testset "string constructor" begin
            # not all symbols twice
            @test_throws ArgumentError arc_diagram(Directed, "AbC,aBd")
            @test_throws ArgumentError arc_diagram(Directed, "ABCD,abdd")
            @test_throws ArgumentError arc_diagram(Directed, "AbdD,aBdD")

            # not all symbols once upper- and once lowercase
            @test_throws ArgumentError arc_diagram(Directed, "ABC,ABC")
            @test_throws ArgumentError arc_diagram(Directed, "abc,aBc")
            @test_throws ArgumentError arc_diagram(Directed, "abDe,BdeA")


            # correct
            @test arc_diagram(Directed, "aACcE,egGIi") ==
                  arc_diagram(Directed, "aACcE", "egGIi") ==
                  arc_diagram(Directed, [1, 0, 0, 1, 0], [0, 0, 1, 1, 0], [-2, -1, -4, -3, 1], [-5, 3, 2, 5, 4])
            @test arc_diagram(Directed, "abcd,ABCD") ==
                  arc_diagram(Directed, "abcd", "ABCD") ==
                  arc_diagram(Directed, [1, 1, 1, 1], [1, 1, 1, 1], [1, 2, 3, 4], [-1, -2, -3, -4])
            @test arc_diagram(Directed, "ABCD,abcd") ==
                  arc_diagram(Directed, "ABCD", "abcd") ==
                  arc_diagram(Directed, [0, 0, 0, 0], [0, 0, 0, 0], [1, 2, 3, 4], [-1, -2, -3, -4])
            @test arc_diagram(Directed, "aBcd,AbCD") ==
                  arc_diagram(Directed, "aBcd", "AbCD") ==
                  arc_diagram(Directed, [1, 0, 1, 1], [1, 0, 1, 1], [1, 2, 3, 4], [-1, -2, -3, -4])
            @test arc_diagram(Directed, "AbBd,aD") ==
                  arc_diagram(Directed, "AbBd", "aD") ==
                  arc_diagram(Directed, [0, 1, 0, 1], [0, 1], [1, -3, -2, 2], [-1, -4])
        end

        @testset "is_crossing_free" begin
            function test_crossing_free(a::ArcDiagramDirected)
                @test is_crossing_free(a) == is_crossing_free(arc_diagram(Undirected, a))
                @test is_crossing_free(a; part=:upper) == is_crossing_free(arc_diagram(Undirected, a); part=:upper)
                @test is_crossing_free(a; part=:lower) == is_crossing_free(arc_diagram(Undirected, a); part=:lower)
            end

            for n in 1:6
                for k in 0:n
                    map(test_crossing_free, all_arc_diagrams(Directed, k, n - k))
                end
            end
        end

        @testset "all_arc_diagrams(Directed, ...)" begin
            @test Set(all_arc_diagrams(Directed, [1, 1, 1, 1], [0, 0])) == Set()
            @test Set(all_arc_diagrams(Directed, [1, 1, 0, 0], [0, 0])) == Set()

            @test Set(all_arc_diagrams(Directed, [1, 0, 0, 1], [1, 0])) == Set(
                arc_diagram(Directed, s) for s in ["aACc,Ee", "aACd,Dc", "aBAb,Ee", "aBAd,Db", "aBCb,Ac", "aBCc,Ab"]
            )

            for n in 0:3
                for i in 0:2n
                    @test length(collect(all_arc_diagrams(Directed, i, 2n - i))) ==
                          length(all_arc_diagrams(Directed, i, 2n - i)) ==
                          2^n * div(factorial(2 * n), 2^n * factorial(n))
                end
                for i in 0:2n+1
                    @test length(collect(all_arc_diagrams(Directed, i, 2n + 1 - i))) ==
                          length(all_arc_diagrams(Directed, i, 2n + 1 - i)) ==
                          0
                end

                for i in 0:2n
                    # two disjoint independent sets of size n each -> matchings on a bipartite graph
                    indep_sets = [[-1:-1:-min(i, n); 1:n-min(i, n)], [-n-1:-1:-max(i, n); n-min(i, n)+1:2n-i]]
                    @test length(collect(all_arc_diagrams(Directed, i, 2n - i; indep_sets))) ==
                          length(all_arc_diagrams(Directed, i, 2n - i; indep_sets)) ==
                          2^n * factorial(n)
                end
            end

            for (parity, num) in [
                ([1, 1, 1, 1, 1, 1], factorial(6)),
                ([0, 1, 1, 1, 1, 1], 3 * 24 * factorial(4) + factorial(6)),
                ([1, 0, 1, 1, 1, 1], 3 * 24 * factorial(4) + factorial(6)),
                ([1, 1, 1, 0, 1, 1], 3 * 24 * factorial(4) + factorial(6)),
                ([0, 0, 1, 1, 1, 1], 6 * binomial(6, 2) * 8 * factorial(2) + 2 * 3 * 24 * factorial(4) + factorial(6)),
                ([0, 1, 1, 0, 1, 1], 4 * binomial(6, 2) * 8 * factorial(2) + 4 * 24 * factorial(4) + factorial(6)),
                ([0, 1, 1, 0, 0, 0], 6 * binomial(6, 2) * 8 * factorial(2) + 2 * 3 * 24 * factorial(4) + factorial(6)),
                ([1, 0, 0, 0, 0, 0], 3 * 24 * factorial(4) + factorial(6)),
                ([0, 0, 0, 0, 0, 0], factorial(6)),
            ]
                iter = all_arc_diagrams(
                    Directed,
                    parity,
                    6,
                    indep_sets=[[-1, -2, -3], [-4, -5, -6], [1, 2], [3, 4], [5, 6]],
                )
                @test length(iter) == length(collect(iter)) == num
            end
        end

    end
end
