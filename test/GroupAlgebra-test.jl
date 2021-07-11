

@testset ExtendedTestSet "All PBWDeformations.GroupAlgebra tests" begin
    @testset "different groupAlgebra constructors" begin
        @testset "not a group" begin
            @test_throws AssertionError PD._groupAlgebra(toGAP(42), "")
            @test_throws AssertionError PD._groupAlgebra(toGAP("Hello World!"), "")
            @test_throws AssertionError PD._groupAlgebra(GAP.SimpleLieAlgebra(toGAP("A"), 2, GAP.Rationals), "")

            @test_throws AssertionError PD._groupAlgebra(GAP.FreeGroup(toGAP("a"), toGAP("b")), "")
            @test_throws AssertionError PD._groupAlgebra(GAP.GL(2, GAP.Integers), "")
        end

        @testset "trivial group" begin
            for ga in [
                    PD.groupAlgebraCyclicGroup(1),
                    PD.groupAlgebraAlternatingGroup(1),
                    PD.groupAlgebraAlternatingGroup(2),
                    PD.groupAlgebraSymmetricGroup(1),
                    ]
                @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == 1
                @test ga.extraData.permRep == ["()"]
                @test ga.relTable == Dict((grp(1), grp(1)) => [(1, [grp(1)])])
            end
        end

        @testset "cyclic group C$n" for n in dimRandomTests
            ga = PD.groupAlgebraCyclicGroup(n)

            @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == n

            showOutput = @test_nowarn sprint(show, ga)
            @test occursin("group algebra", lowercase(showOutput))
            @test occursin("c$n", lowercase(showOutput))

            x = ga.x

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:ga.extraData.order, rand(2:10)))

                # abelian test
                @test PD.normalForm(ga, prod(map(x, ind))) == PD.normalForm(ga, prod(map(x, shuffle(ind))))

                # cyclic test
                @test PD.normalForm(ga, prod(map(x, ind))) == x(1 + sum(map(i -> i-1, ind)) % n)
            end
        end

        @testset "dihedral group D$n" for n in dimRandomTests
            if n % 2 != 0
                @test_throws AssertionError ga = PD.groupAlgebraDihedralGroup(n)
            else
                ga = PD.groupAlgebraDihedralGroup(n)

                @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == n

                showOutput = @test_nowarn sprint(show, ga)
                @test occursin("group algebra", lowercase(showOutput))
                @test occursin("d$n", lowercase(showOutput))

                #TODO: property tests
            end
        end

        @testset "dicyclic group Dic$n" for n in dimRandomTests
            if n % 4 != 0
                @test_throws AssertionError ga = PD.groupAlgebraDicyclicGroup(n)
            else
                ga = PD.groupAlgebraDicyclicGroup(n)

                @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == n

                showOutput = @test_nowarn sprint(show, ga)
                @test occursin("group algebra", lowercase(showOutput))
                @test occursin("dic$n", lowercase(showOutput))

                #TODO: property tests
            end
        end

        @testset "symmetric group S$n" for n in 2:6
            ga = PD.groupAlgebraSymmetricGroup(n)

            @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == factorial(n)

            showOutput = @test_nowarn sprint(show, ga)
            @test occursin("group algebra", lowercase(showOutput))
            @test occursin("s$n", lowercase(showOutput))

            x = ga.x

            # abelian test
            if n <= 2
                for _ in 1:numRandomTests
                    ind = shuffle(rand(1:ga.extraData.order, rand(2:10)))
                    @test PD.normalForm(ga, prod(map(x, ind))) == PD.normalForm(ga, prod(map(x, shuffle(ind))))
                end
            end

            #TODO: more property tests
        end

        @testset "alternating group A$n" for n in 2:6
            ga = PD.groupAlgebraAlternatingGroup(n)

            @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == factorial(n)/2

            showOutput = @test_nowarn sprint(show, ga)
            @test occursin("group algebra", lowercase(showOutput))
            @test occursin("a$n", lowercase(showOutput))

            x = ga.x

            # abelian test
            if n <= 3
                for _ in 1:numRandomTests
                    ind = shuffle(rand(1:ga.extraData.order, rand(2:10)))
                    @test PD.normalForm(ga, prod(map(x, ind))) == PD.normalForm(ga, prod(map(x, shuffle(ind))))
                end
            end

            #TODO: more property tests
        end
    end
    
end
