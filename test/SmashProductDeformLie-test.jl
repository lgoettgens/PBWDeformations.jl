@testset ExtendedTestSet "All SmashProductDeformLie.jl tests" begin
    
    @testset "smash_product_deform_lie constructor" begin
        @testset "$(dynkin)_$n with hw $lambda; R = $R" for (dynkin, n, lambda) in [('A', 2, [1,1]), ('B', 2, [1,0])], R in [QQ, PolynomialRing(QQ, ["x","y","z"])[1]]
            sp, (sp_baseL, sp_baseV) = PD.smash_product_lie(R, dynkin, n, lambda)
            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1,2] = sp_baseL[1]
            kappa[2,1] = -kappa[1,2]
            kappa[3,4] = sp_baseL[2]
            kappa[4,3] = -kappa[3,4]
            deform, (baseL, baseV) = PD.smash_product_deform_lie(sp, kappa)

            @test deform.dimL == sp.dimL == length(deform.baseL)
            @test deform.dimV == sp.dimV == length(deform.baseV)
            @test deform.baseL == baseL
            @test deform.baseV == baseV
            @test deform.coeff_ring == R
            @test deform.symmetric == false
            @test deform.kappa == kappa

            # Test the module basis relations
            for i in 1:length(baseV), j in 1:length(baseV)
                if i == 1 && j == 2
                    @test comm(baseV[i], baseV[j]; strict=true) == baseL[1]
                elseif i == 2 && j == 1
                    @test comm(baseV[i], baseV[j]; strict=true) == -baseL[1]
                elseif i == 3 && j == 4
                    @test comm(baseV[i], baseV[j]; strict=true) == baseL[2]
                elseif i == 4 && j == 3
                    @test comm(baseV[i], baseV[j]; strict=true) == -baseL[2]
                else
                    @test iszero(comm(baseV[i], baseV[j]; strict=true))
                end
            end

            showOutput = @test_nowarn sprint(show, deform)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            # @test_broken occursin(string(dynkin), showOutput)
            # @test_broken occursin(string(lambda), showOutput)
            @test !occursin("symmetric", lowercase(showOutput))
            @test occursin("deformation", lowercase(showOutput))
        end

    end

    @testset "smash_product_symmdeform_lie constructor" begin
        @testset "$(dynkin)_$n with hw $lambda; R = $R" for (dynkin, n, lambda) in [('A', 2, [1,1]), ('B', 2, [1,0])], R in [QQ, PolynomialRing(QQ, ["x","y","z"])[1]]
            sp, _ = PD.smash_product_lie(R, dynkin, n, lambda)
            deform, (baseL, baseV) = PD.smash_product_symmdeform_lie(sp)

            @test deform.dimL == sp.dimL == length(deform.baseL)
            @test deform.dimV == sp.dimV == length(deform.baseV)
            @test deform.baseL == baseL
            @test deform.baseV == baseV
            @test deform.coeff_ring == R
            @test deform.symmetric == true
            @test deform.kappa == fill(zero(sp.alg), sp.dimV, sp.dimV)

            # Test that the module basis commutes
            for vi in baseV, vj in baseV
                @test iszero(comm(vi, vj; strict=true))
            end

            showOutput = @test_nowarn sprint(show, deform)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            # @test_broken occursin(string(dynkin), showOutput)
            # @test_broken occursin(string(lambda), showOutput)
            @test occursin("symmetric deformation", lowercase(showOutput))
        end

    end

    @testset "smash_product_deform_lie sanitize checks; R = $R" for R in [QQ, PolynomialRing(QQ, ["x","y","z"])[1]]
        @testset "check dimensions of kappa" begin
            sp, _ = PD.smash_product_lie(R, 'B', 2, [1,0])

            for eps in [-1,1]
                kappa = fill(zero(sp.alg), sp.dimV+eps, sp.dimV)
                @test_throws ArgumentError("kappa has wrong dimensions.") PD.smash_product_deform_lie(sp, kappa)

                kappa = fill(zero(sp.alg), sp.dimV, sp.dimV+eps)
                @test_throws ArgumentError("kappa has wrong dimensions.") PD.smash_product_deform_lie(sp, kappa)

                kappa = fill(zero(sp.alg), sp.dimV+eps, sp.dimV+eps)
                @test_throws ArgumentError("kappa has wrong dimensions.") PD.smash_product_deform_lie(sp, kappa)
            end
        end

        @testset "check entries of kappa contained in Hopf algebra of smash product" begin
            sp, (baseL, baseV) = PD.smash_product_lie(R, 'B', 2, [1,0])

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1,2] = baseV[1]
            kappa[2,1] = -kappa[1,2]
            @test_throws ArgumentError("kappa does not only take values in the hopf algebra") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1,2] = baseV[1]*baseL[1]
            kappa[2,1] = -kappa[1,2]
            @test_throws ArgumentError("kappa does not only take values in the hopf algebra") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1,2] = sp.alg(2)
            kappa[2,1] = -kappa[1,2]
            @test_nowarn PD.smash_product_deform_lie(sp, kappa)
        end

        @testset "check kappa is skew symmetric" begin
            sp, (baseL, baseV) = PD.smash_product_lie(R, 'B', 2, [1,0])

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1,1] = baseL[1]
            @test_throws ArgumentError("kappa is not skew-symmetric.") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1,2] = baseL[1]
            @test_throws ArgumentError("kappa is not skew-symmetric.") PD.smash_product_deform_lie(sp, kappa)

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1,2] = baseL[1]
            kappa[2,1] = -2*baseL[1]
            @test_throws ArgumentError("kappa is not skew-symmetric.") PD.smash_product_deform_lie(sp, kappa)
        end

    end


    @testset "ispbwdeform" begin
        @testset "symmetric deformation of $(dynkin)_$n with hw $lambda is PBW" for (dynkin, n, lambda) in [('A', 2, [1,1]), ('B', 2, [1,0])]
            sp, _ = PD.smash_product_lie(QQ, dynkin, n, lambda)
            d, _ = PD.smash_product_symmdeform_lie(sp)
            @test PD.ispbwdeform(d)
        end

        @testset "non-PBW deformations" begin
            sp, (baseL, baseV) = PD.smash_product_lie(QQ, 'A', 2, [1,0])
            kappa = fill(zero(sp.alg), 3, 3)
            # some made-up skew-symmetric entries
            kappa[1,2] = baseL[2]
            kappa[2,1] = -baseL[2]
            d, _ = PD.smash_product_deform_lie(sp, kappa)
            @test !PD.ispbwdeform(d)
        end

    end

    @testset "possible_pbwdeforms construction stuff" begin
        @testset "possible_pbwdeforms_nvars tests" begin
            for _ in 1:num_random_tests
                a, b, c = PD.possible_pbwdeforms_nvars(rand(1:10), rand(1:10), rand(0:5))
                @test a == b * c
            end

            for i in 1:10
                @test PD.possible_pbwdeforms_nvars(1, i, 0) == (div(i*(i-1), 2), div(i*(i-1), 2), 1)
                @test PD.possible_pbwdeforms_nvars(1, i, 1) == (2*div(i*(i-1), 2), div(i*(i-1), 2), 2)
                @test PD.possible_pbwdeforms_nvars(5, i, 1) == (6*div(i*(i-1), 2), div(i*(i-1), 2), 6)
                @test PD.possible_pbwdeforms_nvars(1, i, 2) == (3*div(i*(i-1), 2), div(i*(i-1), 2), 3)
                @test PD.possible_pbwdeforms_nvars(3, i, 2) == (10*div(i*(i-1), 2), div(i*(i-1), 2), 10)
            end
        end

        @testset "possible_pbwdeforms_vars tests" begin
            for _ in 1:num_random_tests
                nL, nV, maxdeg = rand(1:10), rand(1:10), rand(0:5)
                l, _, _ = PD.possible_pbwdeforms_nvars(nL, nV, maxdeg)
                @test l == length(PD.possible_pbwdeforms_vars(nL, nV, maxdeg))
            end
        end

        @testset "coefficient_comparison tests" begin
            A, (x,y,z) = PD.free_algebra(QQ, ["x", "y", "z"])
            eq = QQ(2//3)*x + 88*y*z - 12*x*z + 3*y + 0*z^4 - 2*y + 12*x*z
            @test issetequal(PD.coefficient_comparison(eq), elem_type(QQ)[2//3, 88, 1])
        end

        @testset "possible_pbwdeforms tests" begin
            sp, _ = PD.smash_product_lie(QQ, 'A', 1, [1])
            @testset "A_1 with hw [1], maxdeg = $maxdeg" for maxdeg in 0:8
                base = PD.possible_pbwdeforms(sp, maxdeg)

                @test length(base) == 1 + div(maxdeg, 2)

                for b in base
                    for i in 1:sp.dimV, j in 1:sp.dimV
                        @test iszero(b[i,j] + b[j,i])
                    end
                end

                @test repr("text/plain", base[1][1,2]) == "1"
                if length(base) >= 2
                    @test repr("text/plain", base[2][1,2]) == "-2*x_3 + 4*x_1x_2 + x_3^2"
                end
                if length(base) >= 3
                    @test repr("text/plain", base[3][1,2]) == "8*x_3 + 16*x_1x_2 + -32*x_1x_2x_3 + -4*x_3^3 + 16*x_1^2x_2^2 + 8*x_1x_2x_3^2 + x_3^4"
                end
                if length(base) >= 4
                        @test repr("text/plain", base[4][1,2]) == "-96*x_3 + 64*x_1x_2 + -64*x_1x_2x_3 + 40*x_3^3 + 320*x_1^2x_2^2 + 208*x_1x_2x_3^2 + -288*x_1^2x_2^2x_3 + -96*x_1x_2x_3^3 + -6*x_3^5 + 64*x_1^3x_2^3 + 48*x_1^2x_2^2x_3^2 + 12*x_1x_2x_3^4 + x_3^6"
                end
            end

            sp, _ = PD.smash_product_lie(QQ, 'B', 2, [1,0])
            @testset "B_2 with hw [1,0], maxdeg = $maxdeg" for maxdeg in 0:2
                base = PD.possible_pbwdeforms(sp, maxdeg)

                @test length(base) == div(maxdeg+1, 2)

                for b in base
                    for i in 1:sp.dimV, j in 1:sp.dimV
                        @test iszero(b[i,j] + b[j,i])
                    end
                end

                if length(base) >= 1
                    @test repr("text/plain", base[1]) == """
                        5×5 Matrix{PBWDeformations.QuadraticQuoAlgebraElem{fmpq}}:
                         0                    -1*x_4     x_3     -1*x_1      x_9 + 1//2*x_10
                         x_4                  0          x_2     -1//2*x_10  x_5
                         -1*x_3               -1*x_2     0       -1*x_6      x_7
                         x_1                  1//2*x_10  x_6     0           x_8
                         -1*x_9 + -1//2*x_10  -1*x_5     -1*x_7  -1*x_8      0"""
                end
                
            end

        end

    end

       
end
