@testset ExtendedTestSet "All PBWDeformations.SmashProductDeformLie tests" begin
    @testset "smashProducDeformLie coincides with smash_product_symm_deform_lie on symmetric kappa" begin
        @testset "B_2 with hw [1,0]" begin
            sp = PD.smash_product_lie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(AlgebraElement{DefaultScalarType}(), nV, nV)

            deform1 = PD.smash_product_symm_deform_lie(sp)
            deform2 = PD.smash_product_deform_lie(sp, kappa)

            @test deform1 == deform2
        end
    end

    @testset "smash_product_symm_deform_lie constructor" begin
        @testset "$(dynkin)_$n with hw $lambda" for (dynkin, n, lambda) in [('A', 2, [1,1]), ('B', 2, [1,0])]
            sp = PD.smash_product_lie(dynkin, n, lambda)
            deform1 = PD.smash_product_symm_deform_lie(sp)
            deform2 = PD.smash_product_symm_deform_lie(dynkin, n, lambda)

            @test deform1.extraData.symmetric == deform2.extraData.symmetric == true
            @test typeof(deform1.extraData.sp) == typeof(deform2.extraData.sp) == typeof(sp.extraData)
            @test deform1.extraData.sp == deform2.extraData.sp == sp.extraData
            @test deform1.basis == deform2.basis == sp.basis

            # Test that the module basis commutes and the other commutators come from the smash product
            @test deform1.relTable == deform2.relTable == Dict(union(sp.relTable, [(mod(i), mod(j)) => AlgebraElement{DefaultScalarType}([(1//1, mod(j,i))]) for i in 1:sp.extraData.nV for j in 1:i-1]))

            @test sprint(show, deform1) == sprint(show, deform2)

            showOutput = @test_nowarn sprint(show, deform1)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            @test occursin(string(dynkin), showOutput)
            @test occursin(string(lambda), showOutput)
            @test occursin("symmetric deformation", lowercase(showOutput))
        end

    end

    @testset "smash_product_deform_lie assertions" begin
        @testset "assert correct dimensions of kappa" begin
            sp = PD.smash_product_lie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(AlgebraElement{DefaultScalarType}(), nV+1, nV)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(AlgebraElement{DefaultScalarType}(), nV, nV+1)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(AlgebraElement{DefaultScalarType}(), nV+1, nV+1)
            @test_throws AssertionError("size of kappa matches module dimension") PD.smash_product_deform_lie(sp, kappa)
        end

        @testset "assert entries of kappa contained in Hopf algebra of smash product" begin
            sp = PD.smash_product_lie('B', 2, [1,0])
            nV = sp.extraData.nV

            # basis of sp consists of (:mod, 1) to (:mod 5) and (:lie, 1) to (:lie, 10)

            kappa = fill(AlgebraElement{DefaultScalarType}(), nV, nV)
            kappa[1,2] = AlgebraElement{DefaultScalarType}(mod(1))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(AlgebraElement{DefaultScalarType}(), nV, nV)
            kappa[1,2] = AlgebraElement{DefaultScalarType}(lie(0))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(AlgebraElement{DefaultScalarType}(), nV, nV)
            kappa[1,2] = AlgebraElement{DefaultScalarType}(lie(11))
            kappa[2,1] = -kappa[1,2]
            @test_throws AssertionError("kappa only takes values in Hopf algebra") PD.smash_product_deform_lie(sp, kappa)
        end

        @testset "assert kappa is skew symmetric" begin
            sp = PD.smash_product_lie('B', 2, [1,0])
            nV = sp.extraData.nV

            kappa = fill(PD.AlgebraElement{DefaultScalarType}(), nV, nV)
            kappa[1,1] = AlgebraElement{DefaultScalarType}(lie(1))
            @test_throws AssertionError("kappa is skew-symmetric") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(PD.AlgebraElement{DefaultScalarType}(), nV, nV)
            kappa[1,2] = AlgebraElement{DefaultScalarType}(lie(1))
            @test_throws AssertionError("kappa is skew-symmetric") PD.smash_product_deform_lie(sp, kappa)
        end

    end

    @testset "ispbwdeform" begin
        @testset "symmetric deformation of $(dynkin)_$n with hw $lambda" for (dynkin, n, lambda) in [('A', 2, [1,1]), ('B', 2, [1,0])]
            d = PD.smash_product_symm_deform_lie(dynkin, n, lambda)
            @test PD.ispbwdeform(d)
        end

        @testset "non-PBW deformations" begin
            sp = PD.smash_product_lie('A', 2, [1,0])
            kappa = fill(AlgebraElement{DefaultScalarType}(0), 3, 3)
            # some made-up skew-symmetric entries
            kappa[1,2] = 1*lie(3)
            kappa[2,1] = -lie(3)
            d = PD.smash_product_deform_lie(sp, kappa)
            @test !PD.ispbwdeform(d)
        end

    end

    @testset "possible_pbw_deforms construction stuff" begin
        @testset "param_deform_number_vars tests" begin
            for _ in 1:numRandomTests
                a, b, c = PD.param_deform_number_vars(rand(1:10), rand(1:10), rand(0:5))
                @test a == b * c
            end

            for i in 1:10
                @test PD.param_deform_number_vars(1, i, 0) == (div(i*(i-1), 2), div(i*(i-1), 2), 1)
                @test PD.param_deform_number_vars(1, i, 1) == (2*div(i*(i-1), 2), div(i*(i-1), 2), 2)
                @test PD.param_deform_number_vars(5, i, 1) == (6*div(i*(i-1), 2), div(i*(i-1), 2), 6)
                @test PD.param_deform_number_vars(1, i, 2) == (3*div(i*(i-1), 2), div(i*(i-1), 2), 3)
                @test PD.param_deform_number_vars(3, i, 2) == (10*div(i*(i-1), 2), div(i*(i-1), 2), 10)
            end
        end

        @testset "param_deform_vars tests" begin
            for _ in 1:numRandomTests
                nL, nV, maxdeg = rand(1:10), rand(1:10), rand(0:5)
                l, _, _ = PD.param_deform_number_vars(nL, nV, maxdeg)
                @test l == length(PD.param_deform_vars(nL, nV, maxdeg))
            end
        end

        @testset "sort_vars tests" begin
            #TODO
        end

        @testset "coefficient_comparison tests" begin
            eq = 2//3*test(1) + 88*test(3,4) - 12*test(1,5) + 3*test(2) + 0*test(4) - 2*test(2) + 12*test(1,5)
            @test issetequal(PD.coefficient_comparison(eq), [2//3, 88, 1])
        end

        @testset "simplifyGen tests" begin
            #TODO (maybe using Singular)
        end

        @testset "both implementations return the same" begin
            sp = PD.smash_product_lie('A',2,[1,1])
            @test PD.possible_pbw_deforms(sp, 1, use_iterators=false) == PD.possible_pbw_deforms(sp, 1, use_iterators=true)
        end

        @testset "everything still works with use_iterators=$use_iterators with special_return=$special_return" for use_iterators in [false, true], special_return in [Nothing, SparseArrays.SparseMatrixCSC]
            @test_nowarn PD.possible_pbw_deforms(PD.smash_product_lie('A',2,[1,0]), 1; use_iterators, special_return)
            @test_nowarn PD.possible_pbw_deforms(PD.smash_product_lie('A',2,[1,0]), 2; use_iterators, special_return)
            @test_nowarn PD.possible_pbw_deforms(PD.smash_product_lie('A',2,[1,1]), 1; use_iterators, special_return)
            @test_nowarn PD.possible_pbw_deforms(PD.smash_product_lie('B',2,[1,0]), 1; use_iterators, special_return)
        end

    end

end
