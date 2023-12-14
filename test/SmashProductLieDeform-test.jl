@testset "SmashProductLieDeform.jl tests" begin

    @testset "SmashProductLieDeform constructor" begin
        @testset "R = $R" for R in [QQ, cyclotomic_field(4)[1]]
            L = special_orthogonal_lie_algebra(R, 4)
            V = exterior_power(standard_module(L), 2)
            sp = smash_product(L, V)

            kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)))
            kappa[1, 2] = gen(sp, 1, :L).alg_elem
            kappa[2, 1] = -kappa[1, 2]
            kappa[3, 4] = gen(sp, 2, :L).alg_elem
            kappa[4, 3] = -kappa[3, 4]
            d = deform(sp, kappa)

            @test (@inferred coefficient_ring(d)) == R
            @test (@inferred base_lie_algebra(d)) == L
            @test (@inferred base_module(d)) == V

            @test dim(L) == ngens(d, :L)
            @test dim(L) == length(gens(d, :L))
            @test dim(L) == 6

            @test dim(V) == ngens(d, :V)
            @test dim(V) == length(gens(d, :V))
            @test dim(V) == 6

            @test dim(L) + dim(V) == ngens(d)
            @test dim(L) + dim(V) == length(gens(d))
            @test dim(L) + dim(V) == 12

            @test get_attribute(d, :is_symmetric, false) == false
            @test d.kappa == kappa
            @test !iszero(d.kappa)

            # Test the module basis relations
            for i in 1:ngens(d, :V), j in 1:ngens(d, :V)
                if i == 1 && j == 2
                    @test comm(gen(d, i, :V), gen(d, j, :V)) == gen(d, 1, :L)
                elseif i == 2 && j == 1
                    @test comm(gen(d, i, :V), gen(d, j, :V)) == -gen(d, 1, :L)
                elseif i == 3 && j == 4
                    @test comm(gen(d, i, :V), gen(d, j, :V)) == gen(d, 2, :L)
                elseif i == 4 && j == 3
                    @test comm(gen(d, i, :V), gen(d, j, :V)) == -gen(d, 2, :L)
                else
                    @test iszero(comm(gen(d, i, :V), gen(d, j, :V)))
                end
            end

            showOutput = @test_nowarn sprint(show, d)
            @test !occursin("symmetric", lowercase(showOutput))
        end

    end

    @testset "symmetric_deformation constructor" begin
        @testset "R = $R" for R in [QQ, cyclotomic_field(4)[1]]

            for (sp, dimL, dimV) in [begin
                L = special_orthogonal_lie_algebra(R, 4)
                V = exterior_power(standard_module(L), 2)
                sp = smash_product(L, V)
                return (sp, 6, 6)
            end, begin
                L = general_linear_lie_algebra(R, 4)
                V = symmetric_power(standard_module(L), 2)
                sp = smash_product(L, V)
                return (sp, 16, 10)
            end]

                d = symmetric_deformation(sp)

                @test dim(L) == ngens(d, :L)
                @test dim(L) == length(gens(d, :L))
                @test dim(L) == dimL

                @test dim(V) == ngens(d, :V)
                @test dim(V) == length(gens(d, :V))
                @test dim(V) == dimV

                @test dim(L) + dim(V) == ngens(d)
                @test dim(L) + dim(V) == length(gens(d))
                @test dim(L) + dim(V) == dimL + dimV

                @test get_attribute(d, :is_symmetric) == true
                @test iszero(d.kappa)

                # Test that the module basis commutes
                for vi in gens(d, :V), vj in gens(d, :V)
                    @test iszero(comm(vi, vj))
                end

                showOutput = @test_nowarn sprint(show, d)
                @test occursin("symmetric", lowercase(showOutput))
            end
        end
    end

    @testset "SmashProductLieDeform sanitize checks" begin
        @testset "R = $R" for R in [QQ, cyclotomic_field(4)[1]]

            L = special_orthogonal_lie_algebra(R, 4)
            V = exterior_power(standard_module(L), 2)
            sp = smash_product(L, V)

            @testset "check dimensions of kappa" begin
                for eps in [-1, 1]
                    kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)) + eps, dim(base_module(sp)))
                    @test_throws ArgumentError("kappa has wrong dimensions.") deform(sp, kappa)

                    kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)) + eps)
                    @test_throws ArgumentError("kappa has wrong dimensions.") deform(sp, kappa)

                    kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)) + eps, dim(base_module(sp)) + eps)
                    @test_throws ArgumentError("kappa has wrong dimensions.") deform(sp, kappa)

                    kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)) + eps, dim(base_module(sp)) - eps)
                    @test_throws ArgumentError("kappa has wrong dimensions.") deform(sp, kappa)
                end
            end

            @testset "check kappa is skew symmetric" begin
                kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)))
                kappa[1, 1] = gen(sp, 1, :L)
                @test_throws ArgumentError("kappa is not skew-symmetric.") deform(sp, kappa)

                kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)))
                kappa[1, 2] = gen(sp, 1, :L)
                @test_throws ArgumentError("kappa is not skew-symmetric.") deform(sp, kappa)

                kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)))
                kappa[1, 2] = gen(sp, 1, :L)
                kappa[2, 1] = -2 * gen(sp, 1, :L)
                @test_throws ArgumentError("kappa is not skew-symmetric.") deform(sp, kappa)
            end

            @testset "check entries of kappa contained in Hopf algebra of smash product" begin
                kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)))
                kappa[1, 2] = gen(sp, 1, :V)
                kappa[2, 1] = -kappa[1, 2]
                @test_throws ArgumentError("kappa does not only take values in the hopf algebra") deform(sp, kappa)

                kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)))
                kappa[1, 2] = gen(sp, 1, :V) * gen(sp, 1, :L)
                kappa[2, 1] = -kappa[1, 2]
                @test_throws ArgumentError("kappa does not only take values in the hopf algebra") deform(sp, kappa)
            end

            @testset "correct input" begin
                kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)))
                kappa[1, 2] = underlying_algebra(sp)(2)
                kappa[2, 1] = -kappa[1, 2]
                @test_nowarn deform(sp, kappa)

                kappa = zero_matrix(underlying_algebra(sp), dim(base_module(sp)), dim(base_module(sp)))
                kappa[1, 2] = gen(sp, 1, :L) * gen(sp, 2, :L) - 3 * gen(sp, 3, :L)
                kappa[2, 1] = -kappa[1, 2]
                @test_nowarn deform(sp, kappa)
            end
        end
    end

end
