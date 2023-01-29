@testset ExtendedTestSet "All SmashProductDeformLie.jl tests" begin

    @testset "smash_product_deform_lie constructor" begin
        @testset "$(dynkin)_$n with hw $lambda; R = $R" for (dynkin, n, lambda) in [('A', 2, [1, 1]), ('B', 2, [1, 0])],
            R in [QQ, PolynomialRing(QQ, ["x", "y", "z"])[1]]

            sp, (sp_basisL, sp_basisV) = smash_product_lie_highest_weight(R, dynkin, n, lambda)
            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1, 2] = sp_basisL[1]
            kappa[2, 1] = -kappa[1, 2]
            kappa[3, 4] = sp_basisL[2]
            kappa[4, 3] = -kappa[3, 4]
            deform, (basisL, basisV) = smash_product_deform_lie(sp, kappa)

            @test deform.dimL == sp.dimL == length(deform.basisL)
            @test deform.dimV == sp.dimV == length(deform.basisV)
            @test deform.basisL == basisL
            @test deform.basisV == basisV
            @test deform.coeff_ring == R
            @test deform.symmetric == false
            @test deform.kappa == kappa

            @test ngens(deform) == (deform.dimL, deform.dimV)
            @test gens(deform) == (deform.basisL, deform.basisV)

            # Test the module basis relations
            for i in eachindex(basisV), j in eachindex(basisV)
                if i == 1 && j == 2
                    @test normal_form(comm(basisV[i], basisV[j]), deform.rels) == basisL[1]
                elseif i == 2 && j == 1
                    @test normal_form(comm(basisV[i], basisV[j]), deform.rels) == -basisL[1]
                elseif i == 3 && j == 4
                    @test normal_form(comm(basisV[i], basisV[j]), deform.rels) == basisL[2]
                elseif i == 4 && j == 3
                    @test normal_form(comm(basisV[i], basisV[j]), deform.rels) == -basisL[2]
                else
                    @test iszero(normal_form(comm(basisV[i], basisV[j]), deform.rels))
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
        @testset "$(dynkin)_$n with hw $lambda; R = $R" for (dynkin, n, lambda) in [('A', 2, [1, 1]), ('B', 2, [1, 0])],
            R in [QQ, PolynomialRing(QQ, ["x", "y", "z"])[1]]

            sp, _ = smash_product_lie_highest_weight(R, dynkin, n, lambda)
            deform, (basisL, basisV) = smash_product_symmdeform_lie(sp)

            @test deform.dimL == sp.dimL == length(deform.basisL)
            @test deform.dimV == sp.dimV == length(deform.basisV)
            @test deform.basisL == basisL
            @test deform.basisV == basisV
            @test deform.coeff_ring == R
            @test deform.symmetric == true
            @test deform.kappa == fill(zero(sp.alg), sp.dimV, sp.dimV)

            # Test that the module basis commutes
            for vi in basisV, vj in basisV
                @test iszero(normal_form(comm(vi, vj), deform.rels))
            end

            showOutput = @test_nowarn sprint(show, deform)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            # @test_broken occursin(string(dynkin), showOutput)
            # @test_broken occursin(string(lambda), showOutput)
            @test occursin("symmetric deformation", lowercase(showOutput))
        end

    end

    @testset "smash_product_deform_lie sanitize checks; R = $R" for R in [QQ, PolynomialRing(QQ, ["x", "y", "z"])[1]]
        @testset "check dimensions of kappa" begin
            sp, _ = smash_product_lie_highest_weight(R, 'B', 2, [1, 0])

            for eps in [-1, 1]
                kappa = fill(zero(sp.alg), sp.dimV + eps, sp.dimV)
                @test_throws ArgumentError("kappa has wrong dimensions.") smash_product_deform_lie(sp, kappa)

                kappa = fill(zero(sp.alg), sp.dimV, sp.dimV + eps)
                @test_throws ArgumentError("kappa has wrong dimensions.") smash_product_deform_lie(sp, kappa)

                kappa = fill(zero(sp.alg), sp.dimV + eps, sp.dimV + eps)
                @test_throws ArgumentError("kappa has wrong dimensions.") smash_product_deform_lie(sp, kappa)
            end
        end

        @testset "check entries of kappa contained in Hopf algebra of smash product" begin
            sp, (basisL, basisV) = smash_product_lie_highest_weight(R, 'B', 2, [1, 0])

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1, 2] = basisV[1]
            kappa[2, 1] = -kappa[1, 2]
            @test_throws ArgumentError("kappa does not only take values in the hopf algebra") smash_product_deform_lie(
                sp,
                kappa,
            )

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1, 2] = basisV[1] * basisL[1]
            kappa[2, 1] = -kappa[1, 2]
            @test_throws ArgumentError("kappa does not only take values in the hopf algebra") smash_product_deform_lie(
                sp,
                kappa,
            )

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1, 2] = sp.alg(2)
            kappa[2, 1] = -kappa[1, 2]
            @test_nowarn smash_product_deform_lie(sp, kappa)
        end

        @testset "check kappa is skew symmetric" begin
            sp, (basisL, basisV) = smash_product_lie_highest_weight(R, 'B', 2, [1, 0])

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1, 1] = basisL[1]
            @test_throws ArgumentError("kappa is not skew-symmetric.") smash_product_deform_lie(sp, kappa)

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1, 2] = basisL[1]
            @test_throws ArgumentError("kappa is not skew-symmetric.") smash_product_deform_lie(sp, kappa)

            kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
            kappa[1, 2] = basisL[1]
            kappa[2, 1] = -2 * basisL[1]
            @test_throws ArgumentError("kappa is not skew-symmetric.") smash_product_deform_lie(sp, kappa)
        end

    end

end
