C = Rational{Int64}

@testset ExtendedTestSet "All PBWDeformations.GroupAlgebra tests" begin
    @testset "different group_algebra constructors" begin
        @testset "not a group" begin
            @test_throws AssertionError PD.group_algebra(toGAP(42), "")
            @test_throws AssertionError PD.group_algebra(toGAP("Hello World!"), "")
            @test_throws AssertionError PD.group_algebra(GAP.SimpleLieAlgebra(toGAP("A"), 2, GAP.Rationals), "")

            @test_throws AssertionError PD.group_algebra(GAP.FreeGroup(toGAP("a"), toGAP("b")), "")
            @test_throws AssertionError PD.group_algebra(GAP.GL(2, GAP.Integers), "")
        end

        @testset "trivial group" begin
            for ga in [
                    PD.group_algebra_cyclic_group(1),
                    PD.group_algebra_alternating_group(1),
                    PD.group_algebra_alternating_group(2),
                    PD.group_algebra_symmetric_group(1),
                    ]
                @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == 1
                @test ga.extraData.permRep == ["()"]
                @test ga.relTable == Dict((grp(1; C), grp(1; C)) => AlgebraElement{C}([(C(1), grp([1]; C))]))
            end
        end

        @testset "cyclic group C$n" for n in dimRandomTests
            ga = PD.group_algebra_cyclic_group(n)

            @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == n

            showOutput = @test_nowarn sprint(show, ga)
            @test occursin("group algebra", lowercase(showOutput))
            @test occursin("c$n", lowercase(showOutput))

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:ga.extraData.order, rand(2:10)))

                # abelian test
                @test normal_form(ga, grp(ind; C)) ≐ normal_form(ga, grp(shuffle(ind); C))

                # cyclic test
                @test normal_form(ga, grp(ind; C)) ≐ grp(1 + sum(map(i -> i-1, ind)) % n; C)
            end
        end

        @testset "dihedral group D$n" for n in dimRandomTests
            if n % 2 != 0
                @test_throws AssertionError ga = PD.group_algebra_dihedral_group(n)
            else
                ga = PD.group_algebra_dihedral_group(n)

                @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == n

                showOutput = @test_nowarn sprint(show, ga)
                @test occursin("group algebra", lowercase(showOutput))
                @test occursin("d$n", lowercase(showOutput))

                #TODO: property tests
            end
        end

        @testset "dicyclic group Dic$n" for n in dimRandomTests
            if n % 4 != 0
                @test_throws AssertionError ga = PD.group_algebra_dicyclic_group(n)
            else
                ga = PD.group_algebra_dicyclic_group(n)

                @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == n

                showOutput = @test_nowarn sprint(show, ga)
                @test occursin("group algebra", lowercase(showOutput))
                @test occursin("dic$n", lowercase(showOutput))

                #TODO: property tests
            end
        end

        @testset "symmetric group S$n" for n in 2:6
            ga = PD.group_algebra_symmetric_group(n)

            @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == factorial(n)

            showOutput = @test_nowarn sprint(show, ga)
            @test occursin("group algebra", lowercase(showOutput))
            @test occursin("s$n", lowercase(showOutput))

            # abelian test
            if n <= 2
                for _ in 1:numRandomTests
                    ind = shuffle(rand(1:ga.extraData.order, rand(2:10)))
                    @test normal_form(ga, grp(ind; C)) ≐ normal_form(ga, grp(shuffle(ind); C))
                end
            end

            #TODO: more property tests
        end

        @testset "alternating group A$n" for n in 2:6
            ga = PD.group_algebra_alternating_group(n)

            @test ga.extraData.order == length(ga.basis) == length(ga.extraData.permRep) == factorial(n)/2

            showOutput = @test_nowarn sprint(show, ga)
            @test occursin("group algebra", lowercase(showOutput))
            @test occursin("a$n", lowercase(showOutput))

            # abelian test
            if n <= 3
                for _ in 1:numRandomTests
                    ind = shuffle(rand(1:ga.extraData.order, rand(2:10)))
                    @test normal_form(ga, grp(ind; C)) ≐ normal_form(ga, grp(shuffle(ind); C))
                end
            end

            #TODO: more property tests
        end
    end
    
end
