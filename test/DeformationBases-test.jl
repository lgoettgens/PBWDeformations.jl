@testset "DeformationBases/*.jl tests" begin
    @testset "ArcDiagDeformBasis.jl" begin
        @testset "arcdiag_to_deformationmap(:special_orthogonal, :exterior)" begin
            L = special_orthogonal_lie_algebra(QQ, 4, identity_matrix(QQ, 4))
            T = Val(:special_orthogonal)
            V = exterior_power_obj(standard_module(L), 2)
            sp = smash_product(L, V)

            @testset "not all specialisations are zero" begin
                diag = arc_diagram(Undirected, "ABBD,AD")
                dm = PD.arcdiag_to_deformationmap(T, diag, sp)
                @test !iszero(dm)
            end

            @testset "deformation is equivariant, d = $deg" for deg in 1:3
                for diag in PD.pbw_arc_diagrams(T, V, deg)
                    dm = PD.arcdiag_to_deformationmap(T, diag, sp)
                    d = deform(sp, dm)
                    @test all(iszero, pbwdeform_eqs(d, disabled=[:b, :c, :d]))
                end
            end
        end

        @testset "correctness regression" begin
            @testset "SO_4, ⋀²V" begin
                L = special_orthogonal_lie_algebra(QQ, 4, identity_matrix(QQ, 4))
                V = exterior_power_obj(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, ArcDiagDeformBasis{elem_type(sp)}(sp, 0:0); special_return=SMat)[1] ==
                      matrix(QQ, 0, 0, [])
                @test all_pbwdeformations(sp, ArcDiagDeformBasis{elem_type(sp)}(sp, 0:1); special_return=SMat)[1] ==
                      matrix(QQ, 1, 1, [1])

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{MatElem{SmashProductLieElem{QQFieldElem, QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}}:\n [0 x_4 x_5 -x_2 -x_3 0; -x_4 0 x_6 x_1 0 -x_3; -x_5 -x_6 0 0 x_1 x_2; x_2 -x_1 0 0 x_6 -x_5; x_3 0 -x_1 -x_6 0 x_4; 0 x_3 -x_2 x_5 -x_4 0]"
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, 1, 1, [1])
                @test all_pbwdeformations(sp, b) == collect(b)

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:3)
                @test length(collect(b)) == 3
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, [1 0; 0 -3//2; 0 1])
                ms = all_pbwdeformations(sp, b)
                @test length(ms) == 2
                @test ms[1] == collect(b)[1]
                @test 2 * ms[2] == -3 * collect(b)[2] + 2 * collect(b)[3]
            end

            @testset "SO_5, ⋀²V" begin
                L = special_orthogonal_lie_algebra(QQ, 5, identity_matrix(QQ, 5))
                V = exterior_power_obj(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, ArcDiagDeformBasis{elem_type(sp)}(sp, 0:0); special_return=SMat)[1] ==
                      matrix(QQ, 0, 0, [])
                @test all_pbwdeformations(sp, ArcDiagDeformBasis{elem_type(sp)}(sp, 0:1); special_return=SMat)[1] ==
                      matrix(QQ, 1, 1, [1])

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{MatElem{SmashProductLieElem{QQFieldElem, QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}}:\n [0 x_5 x_6 x_7 -x_2 -x_3 -x_4 0 0 0; -x_5 0 x_8 x_9 x_1 0 0 -x_3 -x_4 0; -x_6 -x_8 0 x_10 0 x_1 0 x_2 0 -x_4; -x_7 -x_9 -x_10 0 0 0 x_1 0 x_2 x_3; x_2 -x_1 0 0 0 x_8 x_9 -x_6 -x_7 0; x_3 0 -x_1 0 -x_8 0 x_10 x_5 0 -x_7; x_4 0 0 -x_1 -x_9 -x_10 0 0 x_5 x_6; 0 x_3 -x_2 0 x_6 -x_5 0 0 x_10 -x_9; 0 x_4 0 -x_2 x_7 0 -x_5 -x_10 0 x_8; 0 0 x_4 -x_3 0 x_7 -x_6 x_9 -x_8 0]"
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, 1, 1, [1])
                @test all_pbwdeformations(sp, b) == collect(b)
            end

            @testset "SO_4, S²V" begin
                L = special_orthogonal_lie_algebra(QQ, 4, identity_matrix(QQ, 4))
                V = symmetric_power_obj(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, ArcDiagDeformBasis{elem_type(sp)}(sp, 0:0); special_return=SMat)[1] ==
                      matrix(QQ, 0, 0, [])
                @test all_pbwdeformations(sp, ArcDiagDeformBasis{elem_type(sp)}(sp, 0:1); special_return=SMat)[1] ==
                      matrix(QQ, 1, 1, [1])

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:2)
                @test length(collect(b)) == 2
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, 2, 1, [1, 0])
                @test all_pbwdeformations(sp, b) == collect(b)[1:1]
            end

            @testset "SO_5, S²V" begin
                L = special_orthogonal_lie_algebra(QQ, 5, identity_matrix(QQ, 5))
                V = symmetric_power_obj(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, ArcDiagDeformBasis{elem_type(sp)}(sp, 0:0); special_return=SMat)[1] ==
                      matrix(QQ, 0, 0, [])
                @test all_pbwdeformations(sp, ArcDiagDeformBasis{elem_type(sp)}(sp, 0:1); special_return=SMat)[1] ==
                      matrix(QQ, 1, 1, [1])

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:2)
                @test length(collect(b)) == 2
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, 2, 1, [1, 0])
                @test all_pbwdeformations(sp, b) == collect(b)[1:1]
            end

            @testset "SO_2, T²V" begin
                L = special_orthogonal_lie_algebra(QQ, 2, identity_matrix(QQ, 2))
                V = tensor_power_obj(standard_module(L), 2)
                sp = smash_product(L, V)

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:0)
                @test length(b) == 0
                ms = all_pbwdeformations(sp, b)
                @test length(ms) == 0
            end

            @testset "SO_3, T²V" begin
                L = special_orthogonal_lie_algebra(QQ, 3, identity_matrix(QQ, 3))
                V = tensor_power_obj(standard_module(L), 2)
                sp = smash_product(L, V)

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:0)
                @test length(b) == 0
                ms = all_pbwdeformations(sp, b)
                @test length(ms) == 0
            end

            @testset "SO_2, T³V" begin
                deformmap(sp::SmashProductLie, diag::String) = PBWDeformations.normalize_default(
                    PBWDeformations.arcdiag_to_deformationmap(PBWDeformations.SO(), arc_diagram(Undirected, diag), sp),
                )

                L = special_orthogonal_lie_algebra(QQ, 2, identity_matrix(QQ, 2))
                V = tensor_power_obj(standard_module(L), 3)
                sp = smash_product(L, V)

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:0)
                @test length(b) == 3
                ms = all_pbwdeformations(sp, b)
                @test length(ms) == 3
                @test issetequal(ms, [
                    deformmap(sp, "AACCEE,"), # same as "ABBDDA,"
                    deformmap(sp, "AACDCD,"), # same as "ABADDB,"
                    deformmap(sp, "ABABEE,"), # same as "ABBDAD,"
                    # deformmap(sp, "ABCCAB,"), # same as "ABCBCA,"; already linear depedent on the others
                ])
            end

            @testset "SO_3, T³V" begin
                deformmap(sp::SmashProductLie, diag::String) = PBWDeformations.normalize_default(
                    PBWDeformations.arcdiag_to_deformationmap(PBWDeformations.SO(), arc_diagram(Undirected, diag), sp),
                )

                L = special_orthogonal_lie_algebra(QQ, 3, identity_matrix(QQ, 3))
                V = tensor_power_obj(standard_module(L), 3)
                sp = smash_product(L, V)

                b = ArcDiagDeformBasis{elem_type(sp)}(sp, 0:0)
                @test length(b) == 4
                ms = all_pbwdeformations(sp, b)
                @test length(ms) == 4
                @test issetequal(
                    ms,
                    [
                        deformmap(sp, "AACCEE,"), # same as "ABBDDA,"
                        deformmap(sp, "AACDCD,"), # same as "ABADDB,"
                        deformmap(sp, "ABABEE,"), # same as "ABBDAD,"
                        deformmap(sp, "ABCCAB,"), # same as "ABCBCA,"
                    ],
                )
            end
        end

    end

    @testset "PseudographDeformBasis.jl" begin
        @testset "correctness regression" begin
            @testset "SO_4, ⋀²V" begin
                L = special_orthogonal_lie_algebra(QQ, 4, identity_matrix(QQ, 4))
                V = exterior_power_obj(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, PseudographDeformBasis{elem_type(sp)}(sp, 0:0); special_return=SMat)[1] ==
                      matrix(QQ, 0, 0, [])
                @test all_pbwdeformations(sp, PseudographDeformBasis{elem_type(sp)}(sp, 0:1); special_return=SMat)[1] ==
                      matrix(QQ, 1, 1, [1])

                b = PseudographDeformBasis{elem_type(sp)}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{MatElem{SmashProductLieElem{QQFieldElem, QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}}:\n [0 x_4 x_5 -x_2 -x_3 0; -x_4 0 x_6 x_1 0 -x_3; -x_5 -x_6 0 0 x_1 x_2; x_2 -x_1 0 0 x_6 -x_5; x_3 0 -x_1 -x_6 0 x_4; 0 x_3 -x_2 x_5 -x_4 0]"
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, 1, 1, [1])
                @test all_pbwdeformations(sp, b) == collect(b)

                b = PseudographDeformBasis{elem_type(sp)}(sp, 0:3)
                @test length(collect(b)) == 3
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, [1 0; 0 -3//2; 0 1])
                ms = all_pbwdeformations(sp, b)
                @test length(ms) == 2
                @test ms[1] == collect(b)[1]
                @test 2 * ms[2] == -3 * collect(b)[2] + 2 * collect(b)[3]
            end

            @testset "SO_5, ⋀²V" begin
                L = special_orthogonal_lie_algebra(QQ, 5, identity_matrix(QQ, 5))
                V = exterior_power_obj(standard_module(L), 2)
                sp = smash_product(L, V)

                @test all_pbwdeformations(sp, PseudographDeformBasis{elem_type(sp)}(sp, 0:0); special_return=SMat)[1] ==
                      matrix(QQ, 0, 0, [])
                @test all_pbwdeformations(sp, PseudographDeformBasis{elem_type(sp)}(sp, 0:1); special_return=SMat)[1] ==
                      matrix(QQ, 1, 1, [1])

                b = PseudographDeformBasis{elem_type(sp)}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{MatElem{SmashProductLieElem{QQFieldElem, QQFieldElem, LinearLieAlgebraElem{QQFieldElem}}}}:\n [0 x_5 x_6 x_7 -x_2 -x_3 -x_4 0 0 0; -x_5 0 x_8 x_9 x_1 0 0 -x_3 -x_4 0; -x_6 -x_8 0 x_10 0 x_1 0 x_2 0 -x_4; -x_7 -x_9 -x_10 0 0 0 x_1 0 x_2 x_3; x_2 -x_1 0 0 0 x_8 x_9 -x_6 -x_7 0; x_3 0 -x_1 0 -x_8 0 x_10 x_5 0 -x_7; x_4 0 0 -x_1 -x_9 -x_10 0 0 x_5 x_6; 0 x_3 -x_2 0 x_6 -x_5 0 0 x_10 -x_9; 0 x_4 0 -x_2 x_7 0 -x_5 -x_10 0 x_8; 0 0 x_4 -x_3 0 x_7 -x_6 x_9 -x_8 0]"
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, 1, 1, [1])
                @test all_pbwdeformations(sp, b) == collect(b)

                b = PseudographDeformBasis{elem_type(sp)}(sp, 0:3)
                @test length(collect(b)) == 4
                @test all_pbwdeformations(sp, b; special_return=SMat)[1] == matrix(QQ, 4, 1, [1, 0, 0, 0])
                @test all_pbwdeformations(sp, b) == collect(b)[1:1]
            end
        end

    end

end
