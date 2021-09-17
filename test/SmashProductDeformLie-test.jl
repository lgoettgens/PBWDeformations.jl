@testset ExtendedTestSet "All PBWDeformations.SmashProductDeformLie tests" begin
    @testset "smashProducDeformLie coincides with smashProductSymmDeformLie on symmetric kappa" begin
        @testset "B_2 with hw [1,0]" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(AlgebraElement{Rational{Int64}}(), nV, nV)

            deform1 = PD.smashProductSymmDeformLie(sp)
            deform2 = PD.smashProductDeformLie(sp, kappa)

            @test deform1 == deform2
        end
    end

    @testset "smashProductSymmDeformLie constructor" begin
        @testset "$(dynkin)_$n with hw $lambda" for (dynkin, n, lambda) in [('A', 2, [1,1]), ('B', 2, [1,0])]
            sp = PD.smashProductLie(dynkin, n, lambda)
            deform1 = PD.smashProductSymmDeformLie(sp)
            deform2 = PD.smashProductSymmDeformLie(dynkin, n, lambda)

            @test deform1.extraData.symmetric == deform2.extraData.symmetric == true
            @test typeof(deform1.extraData.sp) == typeof(deform2.extraData.sp) == typeof(sp.extraData)
            @test deform1.extraData.sp == deform2.extraData.sp == sp.extraData
            @test deform1.basis == deform2.basis == sp.basis

            # Test that the module basis commutes and the other commutators come from the smash product
            @test deform1.relTable == deform2.relTable == Dict(union(sp.relTable, [(mod(i), mod(j)) => AlgebraElement{Rational{Int64}}([(1//1, mod(j,i))]) for i in 1:sp.extraData.nV for j in 1:i-1]))

            @test sprint(show, deform1) == sprint(show, deform2)

            showOutput = @test_nowarn sprint(show, deform1)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            @test occursin(string(dynkin), showOutput)
            @test occursin(string(lambda), showOutput)
            @test occursin("symmetric deformation", lowercase(showOutput))
        end

    end

    @testset "smashProductDeformLie assertions" begin
        @testset "assert correct dimensions of kappa" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(AlgebraElement{Rational{Int64}}(), nV+1, nV)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(AlgebraElement{Rational{Int64}}(), nV, nV+1)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(AlgebraElement{Rational{Int64}}(), nV+1, nV+1)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smashProductDeformLie(sp, kappa)
        end

        @testset "assert entries of kappa contained in Hopf algebra of smash product" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            # basis of sp consists of (:mod, 1) to (:mod 5) and (:lie, 1) to (:lie, 10)

            kappa = fill(AlgebraElement{Rational{Int64}}(), nV, nV)
            kappa[1,2] = AlgebraElement{Rational{Int64}}(mod(1))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(AlgebraElement{Rational{Int64}}(), nV, nV)
            kappa[1,2] = AlgebraElement{Rational{Int64}}(lie(0))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(AlgebraElement{Rational{Int64}}(), nV, nV)
            kappa[1,2] = AlgebraElement{Rational{Int64}}(lie(11))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smashProductDeformLie(sp, kappa)
        end

        @testset "assert kappa is skew symmetric" begin
            sp = PD.smashProductLie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(PD.AlgebraElement{Rational{Int64}}(), nV, nV)
            kappa[1,1] = AlgebraElement{Rational{Int64}}(lie(1))
            @test_throws AssertionError("kappa is skew-symmetric") PD.smashProductDeformLie(sp, kappa)

            kappa = fill(PD.AlgebraElement{Rational{Int64}}(), nV, nV)
            kappa[1,2] = AlgebraElement{Rational{Int64}}(lie(1))
            @test_throws AssertionError("kappa is skew-symmetric") PD.smashProductDeformLie(sp, kappa)
        end

    end

    @testset "isPBWDeform" begin
        @testset "symmetric deformation of $(dynkin)_$n with hw $lambda" for (dynkin, n, lambda) in [('A', 2, [1,1]), ('B', 2, [1,0])]
            d = PD.smashProductSymmDeformLie(dynkin, n, lambda)
            @test PD.isPBWDeform(d)
        end

        @testset "non-PBW deformations" begin
            sp = PD.smashProductLie('A', 2, [1,0])
            kappa = fill(AlgebraElement{Rational{Int64}}(0), 3, 3)
            # some made-up skew-symmetric entries
            kappa[1,2] = 1*lie(3)
            kappa[2,1] = -lie(3)
            d = PD.smashProductDeformLie(sp, kappa)
            @test !PD.isPBWDeform(d)
        end

    end

    @testset "possible_pbw_deforms construction stuff" begin
        @testset "paramDeformNumberVars tests" begin
            for _ in 1:numRandomTests
                a, b, c = PD.paramDeformNumberVars(rand(1:10), rand(1:10), rand(0:5))
                @test a == b * c
            end

            for i in 1:10
                @test PD.paramDeformNumberVars(1, i, 0) == (div(i*(i-1), 2), div(i*(i-1), 2), 1)
                @test PD.paramDeformNumberVars(1, i, 1) == (2*div(i*(i-1), 2), div(i*(i-1), 2), 2)
                @test PD.paramDeformNumberVars(5, i, 1) == (6*div(i*(i-1), 2), div(i*(i-1), 2), 6)
                @test PD.paramDeformNumberVars(1, i, 2) == (3*div(i*(i-1), 2), div(i*(i-1), 2), 3)
                @test PD.paramDeformNumberVars(3, i, 2) == (10*div(i*(i-1), 2), div(i*(i-1), 2), 10)
            end
        end

        @testset "paramDeformVars tests" begin
            for _ in 1:numRandomTests
                nL, nV, maxdeg = rand(1:10), rand(1:10), rand(0:5)
                l, _, _ = PD.paramDeformNumberVars(nL, nV, maxdeg)
                @test l == length(PD.paramDeformVars(nL, nV, maxdeg))
            end
        end

        @testset "sortVars tests" begin
            #TODO
        end

        @testset "coefficientComparison tests" begin
            eq = 2//3*test(1) + 88*test(3,4) - 12*test(1,5) + 3*test(2) + 0*test(4) - 2*test(2) + 12*test(1,5)
            @test issetequal(PD.coefficientComparison(eq), [2//3, 88, 1])
        end

        @testset "simplifyGen tests" begin
            #TODO (maybe using Singular)
        end

        @testset "both implementations return the same" begin
            sp = PD.smashProductLie('A',2,[1,1])
            @test PD.possible_pbw_deforms(sp, 1, use_iterators=false) == PD.possible_pbw_deforms(sp, 1, use_iterators=true)
        end

        @testset "everything still works with use_iterators=$use_iterators with special_return=$special_return" for use_iterators in [false, true], special_return in [Nothing, SparseArrays.SparseMatrixCSC]
            @test_nowarn PD.possible_pbw_deforms(PD.smashProductLie('A',2,[1,0]), 1; use_iterators, special_return)
            @test_nowarn PD.possible_pbw_deforms(PD.smashProductLie('A',2,[1,0]), 2; use_iterators, special_return)
            @test_nowarn PD.possible_pbw_deforms(PD.smashProductLie('A',2,[1,1]), 1; use_iterators, special_return)
            @test_nowarn PD.possible_pbw_deforms(PD.smashProductLie('B',2,[1,0]), 1; use_iterators, special_return)
        end

    end

end
