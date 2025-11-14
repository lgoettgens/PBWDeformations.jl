@testset "SmashProductPBWDeformLie.jl tests" begin
    @testset "is_pbwdeformation" begin
        @testset "symmetric deformation of so_4(QQ) ⋉ ⋀^2 V" begin
            L = special_orthogonal_lie_algebra(QQ, 4, identity_matrix(QQ, 4))
            V = exterior_power_obj(standard_module(L), 2)
            sp = smash_product(L, V)

            d = symmetric_deformation(sp)
            @test is_pbwdeformation(d)
        end

        @testset "non-PBW deformations" begin
            L = special_orthogonal_lie_algebra(QQ, 4, identity_matrix(QQ, 4))
            V = exterior_power_obj(standard_module(L), 2)
            sp = smash_product(L, V)

            kappa = zero_matrix(sp, 6, 6)
            # some made-up skew-symmetric entries
            kappa[1, 2] = gen(sp, :L, 2)
            kappa[2, 1] = -gen(sp, :L, 2)
            d = deform(sp, kappa)
            @test !is_pbwdeformation(d)
        end

    end

    @testset "all_pbwdeformations construction stuff" begin
        @testset "coefficient_comparison tests" begin
            A, (x, y, z) = free_associative_algebra(QQ, ["x", "y", "z"])
            eq = QQ(2 // 3) * x + 88 * y * z - 12 * x * z + 3 * y + 0 * z^4 - 2 * y + 12 * x * z
            @test issetequal(PD.coefficient_comparison(eq), elem_type(QQ)[2//3, 88, 1])
        end

        @testset "all_pbwdeformations tests" begin
            L = lie_algebra(QQ, :A, 1)
            V = simple_module(L, [1])
            sp = smash_product(L, V)
            @testset "A_1 with hw [1], maxdeg = $maxdeg" for maxdeg in 0:8
                basis = all_pbwdeformations(sp, 0:maxdeg)

                @test length(basis) == 1 + div(maxdeg, 2)

                for b in basis
                    for i in 1:dim(base_module(sp)), j in 1:dim(base_module(sp))
                        @test iszero(b[i, j] + b[j, i])
                    end
                end

                @test repr("text/plain", basis[1][1, 2]) == "1"
                if length(basis) >= 2
                    @test repr("text/plain", basis[2][1, 2]) == "4*x_1*y_1 + h_1^2 - 2*h_1"
                end
                if length(basis) >= 3
                    @test repr("text/plain", basis[3][1, 2]) ==
                          "16*x_1^2*y_1^2 + 8*x_1*y_1*h_1^2 + h_1^4 - 32*x_1*y_1*h_1 - 4*h_1^3 + 16*x_1*y_1 + 8*h_1"
                end
                if length(basis) >= 4
                    @test repr("text/plain", basis[4][1, 2]) ==
                          "64*x_1^3*y_1^3 + 48*x_1^2*y_1^2*h_1^2 + 12*x_1*y_1*h_1^4 + h_1^6 - 288*x_1^2*y_1^2*h_1 - 96*x_1*y_1*h_1^3 - 6*h_1^5 + 320*x_1^2*y_1^2 + 208*x_1*y_1*h_1^2 - 64*x_1*y_1*h_1 + 40*h_1^3 + 64*x_1*y_1 - 96*h_1"
                end
            end
        end
    end

end
