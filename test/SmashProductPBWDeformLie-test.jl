normalize_basis = PD.normalize_basis

@testset ExtendedTestSet "All SmashProductPBWDeformLie.jl tests" begin

    @testset "is_pbwdeform" begin
        @testset "symmetric deformation of $(dynkin)_$n with hw $lambda is PBW" for (dynkin, n, lambda) in
                                                                                    [('A', 2, [1, 1]), ('B', 2, [1, 0])]
            sp, _ = smash_product_lie_highest_weight(QQ, dynkin, n, lambda)
            d, _ = smash_product_symmdeform_lie(sp)
            @test is_pbwdeform(d)
        end

        @testset "non-PBW deformations" begin
            sp, (basisL, basisV) = smash_product_lie_highest_weight(QQ, 'A', 2, [1, 0])
            kappa = fill(zero(sp.alg), 3, 3)
            # some made-up skew-symmetric entries
            kappa[1, 2] = basisL[2]
            kappa[2, 1] = -basisL[2]
            d, _ = smash_product_deform_lie(sp, kappa)
            @test !is_pbwdeform(d)
        end

    end

    @testset "pbwdeforms_all construction stuff" begin
        @testset "coefficient_comparison tests" begin
            A, (x, y, z) = free_algebra(QQ, ["x", "y", "z"])
            eq = QQ(2 // 3) * x + 88 * y * z - 12 * x * z + 3 * y + 0 * z^4 - 2 * y + 12 * x * z
            @test issetequal(PD.coefficient_comparison(eq), elem_type(QQ)[2//3, 88, 1])
        end

        @testset "normalize_basis tests" begin
            A, (x, y) = PolynomialRing(QQ, ["x", "y"])
            b0 = [A(0) A(0); A(0) A(0)]
            b1 = [A(0) A(1); A(-1) A(0)]
            b2 = [A(0) x; -x A(0)]
            b3 = [A(0) y; -y A(0)]

            @test issetequal(normalize_basis([b0]), [])
            @test length(normalize_basis([b0])) == 0

            @test issetequal(normalize_basis([b1]), [b1])
            @test length(normalize_basis([b1])) == 1

            @test issetequal(normalize_basis([b2]), [b2])
            @test length(normalize_basis([b2])) == 1

            @test issetequal(normalize_basis([b1, b3]), [b1, b3])
            @test length(normalize_basis([b1, b3])) == 2

            @test issetequal(normalize_basis([b1, b0]), [b1])
            @test length(normalize_basis([b1, b0])) == 1

            @test issetequal(normalize_basis([b1, b1]), [b1])
            @test length(normalize_basis([b1, b1])) == 1

            @test issetequal(normalize_basis([b1, 2 * b1]), [b1])
            @test length(normalize_basis([b1, 2 * b1])) == 1
        end


        @testset "pbwdeforms_all tests" begin
            sp, _ = smash_product_lie_highest_weight(QQ, 'A', 1, [1])
            @testset "A_1 with hw [1], maxdeg = $maxdeg" for maxdeg in 0:8
                basis = pbwdeforms_all(sp, 0:maxdeg)

                @test length(basis) == 1 + div(maxdeg, 2)

                for b in basis
                    for i in 1:sp.dimV, j in 1:sp.dimV
                        @test iszero(b[i, j] + b[j, i])
                    end
                end

                @test repr("text/plain", basis[1][1, 2]) == "1"
                if length(basis) >= 2
                    @test repr("text/plain", basis[2][1, 2]) == "-2*x_3 + 4*x_1*x_2 + x_3^2"
                end
                if length(basis) >= 3
                    @test repr("text/plain", basis[3][1, 2]) ==
                          "8*x_3 + 16*x_1*x_2 + -32*x_1*x_2*x_3 + -4*x_3^3 + 16*x_1^2*x_2^2 + 8*x_1*x_2*x_3^2 + x_3^4"
                end
                if length(basis) >= 4
                    @test repr("text/plain", basis[4][1, 2]) ==
                          "-96*x_3 + 64*x_1*x_2 + -64*x_1*x_2*x_3 + 40*x_3^3 + 320*x_1^2*x_2^2 + 208*x_1*x_2*x_3^2 + -288*x_1^2*x_2^2*x_3 + -96*x_1*x_2*x_3^3 + -6*x_3^5 + 64*x_1^3*x_2^3 + 48*x_1^2*x_2^2*x_3^2 + 12*x_1*x_2*x_3^4 + x_3^6"
                end
            end

            sp, _ = smash_product_lie_highest_weight(QQ, 'B', 2, [1, 0])
            @testset "B_2 with hw [1,0], maxdeg = $maxdeg" for maxdeg in 0:2
                basis = pbwdeforms_all(sp, 0:maxdeg)

                @test length(basis) == div(maxdeg + 1, 2)

                for b in basis
                    for i in 1:sp.dimV, j in 1:sp.dimV
                        @test iszero(b[i, j] + b[j, i])
                    end
                end

                if length(basis) >= 1
                    @test repr("text/plain", basis[1]) == """
                        5×5 Matrix{QuadraticQuoAlgebraElem{fmpq}}:
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
