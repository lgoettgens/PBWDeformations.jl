@testset ExtendedTestSet "All DeformationBases/*.jl tests" begin
    @testset "All DeformationBases/ArcDiagDeformBasis.jl tests" begin
        @testset "arcdiag_to_basiselem__so_powers_stdmod(:exterior)" begin
            sp, _ = smash_product_lie_so_extpowers_standard_module(QQ, 4, 2)
            @testset "not all specialisations are zero" begin
                diag = ArcDiagram(4, 2, [5, 3, 2, 6, 1, 4])
                dm = PD.arcdiag_to_basiselem__so_powers_stdmod(diag, 4, :exterior, 2, 1, sp.alg(0), sp.rels)
                @test !iszero(dm)
            end

            @testset "deformation is equivariant, d = $d" for d in 1:3
                for diag in PD.pbw_arc_diagrams__so_powers_stdmod(:exterior, 2, d)
                    dm = PD.arcdiag_to_basiselem__so_powers_stdmod(diag, 4, :exterior, 2, d, sp.alg(0), sp.rels)
                    deform, _ = smash_product_deform_lie(sp, dm)
                    @test all(iszero, pbwdeform_eqs(deform, disabled=[:b, :c, :d]))
                end
            end
        end

        @testset "correctness regression" begin
            @testset "SO_4, ⋀²V" begin
                sp, _ = smash_product_lie_so_extpowers_standard_module(QQ, 4, 2)

                @test pbwdeforms_all(sp, ArcDiagDeformBasis{fmpq}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 0, 0)
                @test pbwdeforms_all(sp, ArcDiagDeformBasis{fmpq}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = ArcDiagDeformBasis{fmpq}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{Any}:\n AbstractAlgebra.Generic.FreeAssAlgElem{fmpq}[0 x_2_3 x_2_4 -x_1_3 -x_1_4 0; -x_2_3 0 x_3_4 x_1_2 0 -x_1_4; -x_2_4 -x_3_4 0 0 x_1_2 x_1_3; x_1_3 -x_1_2 0 0 x_3_4 -x_2_4; x_1_4 0 -x_1_2 -x_3_4 0 x_2_3; 0 x_1_4 -x_1_3 x_2_4 -x_2_3 0]"
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] == sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = ArcDiagDeformBasis{fmpq}(sp, 0:3)
                @test length(collect(b)) == 4
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse([2, 2, 2], [2, 3, 4], fmpq[1, 3//2, -1//2], 4, 4)
            end

            @testset "SO_5, ⋀²V" begin
                sp, _ = smash_product_lie_so_extpowers_standard_module(QQ, 5, 2)

                @test pbwdeforms_all(sp, ArcDiagDeformBasis{fmpq}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 0, 0)
                @test pbwdeforms_all(sp, ArcDiagDeformBasis{fmpq}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = ArcDiagDeformBasis{fmpq}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{Any}:\n AbstractAlgebra.Generic.FreeAssAlgElem{fmpq}[0 x_2_3 x_2_4 x_2_5 -x_1_3 -x_1_4 -x_1_5 0 0 0; -x_2_3 0 x_3_4 x_3_5 x_1_2 0 0 -x_1_4 -x_1_5 0; -x_2_4 -x_3_4 0 x_4_5 0 x_1_2 0 x_1_3 0 -x_1_5; -x_2_5 -x_3_5 -x_4_5 0 0 0 x_1_2 0 x_1_3 x_1_4; x_1_3 -x_1_2 0 0 0 x_3_4 x_3_5 -x_2_4 -x_2_5 0; x_1_4 0 -x_1_2 0 -x_3_4 0 x_4_5 x_2_3 0 -x_2_5; x_1_5 0 0 -x_1_2 -x_3_5 -x_4_5 0 0 x_2_3 x_2_4; 0 x_1_4 -x_1_3 0 x_2_4 -x_2_3 0 0 x_4_5 -x_3_5; 0 x_1_5 0 -x_1_3 x_2_5 0 -x_2_3 -x_4_5 0 x_3_4; 0 0 x_1_5 -x_1_4 0 x_2_5 -x_2_4 x_3_5 -x_3_4 0]"
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] == sparse(Int64[], Int64[], fmpq[], 1, 1)
            end

            @testset "SO_4, S²V" begin
                sp, _ = smash_product_lie_so_symmpowers_standard_module(QQ, 4, 2)

                @test pbwdeforms_all(sp, ArcDiagDeformBasis{fmpq}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 0, 0)
                @test pbwdeforms_all(sp, ArcDiagDeformBasis{fmpq}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = ArcDiagDeformBasis{fmpq}(sp, 0:2)
                @test length(collect(b)) == 2
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[2], Int64[2], fmpq[1], 2, 2)
            end

            @testset "SO_5, S²V" begin
                sp, _ = smash_product_lie_so_symmpowers_standard_module(QQ, 5, 2)

                @test pbwdeforms_all(sp, ArcDiagDeformBasis{fmpq}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 0, 0)
                @test pbwdeforms_all(sp, ArcDiagDeformBasis{fmpq}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = ArcDiagDeformBasis{fmpq}(sp, 0:2)
                @test length(collect(b)) == 2
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[2], Int64[2], fmpq[1], 2, 2)
            end
        end

    end

    @testset "All DeformationBases/PseudographDeformBasis.jl tests" begin
        @testset "correctness regression" begin
            @testset "SO_4, ⋀²V" begin
                sp, _ = smash_product_lie_so_extpowers_standard_module(QQ, 4, 2)

                @test pbwdeforms_all(sp, PseudographDeformBasis{fmpq}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 0, 0)
                @test pbwdeforms_all(sp, PseudographDeformBasis{fmpq}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = PseudographDeformBasis{fmpq}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{Any}:\n AbstractAlgebra.Generic.FreeAssAlgElem{fmpq}[0 x_2_3 x_2_4 -x_1_3 -x_1_4 0; -x_2_3 0 x_3_4 x_1_2 0 -x_1_4; -x_2_4 -x_3_4 0 0 x_1_2 x_1_3; x_1_3 -x_1_2 0 0 x_3_4 -x_2_4; x_1_4 0 -x_1_2 -x_3_4 0 x_2_3; 0 x_1_4 -x_1_3 x_2_4 -x_2_3 0]"
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] == sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = PseudographDeformBasis{fmpq}(sp, 0:3)
                @test length(collect(b)) == 4
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse([2, 2, 2], [2, 3, 4], fmpq[1, 3//2, -1//2], 4, 4)
            end

            @testset "SO_5, ⋀²V" begin
                sp, _ = smash_product_lie_so_extpowers_standard_module(QQ, 5, 2)

                @test pbwdeforms_all(sp, PseudographDeformBasis{fmpq}(sp, 0:0); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 0, 0)
                @test pbwdeforms_all(sp, PseudographDeformBasis{fmpq}(sp, 0:1); special_return=SparseMatrixCSC)[1] ==
                      sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = PseudographDeformBasis{fmpq}(sp, 0:2)
                @test length(collect(b)) == 1
                @test repr("text/plain", collect(b)) ==
                      "1-element Vector{Any}:\n AbstractAlgebra.Generic.FreeAssAlgElem{fmpq}[0 x_2_3 x_2_4 x_2_5 -x_1_3 -x_1_4 -x_1_5 0 0 0; -x_2_3 0 x_3_4 x_3_5 x_1_2 0 0 -x_1_4 -x_1_5 0; -x_2_4 -x_3_4 0 x_4_5 0 x_1_2 0 x_1_3 0 -x_1_5; -x_2_5 -x_3_5 -x_4_5 0 0 0 x_1_2 0 x_1_3 x_1_4; x_1_3 -x_1_2 0 0 0 x_3_4 x_3_5 -x_2_4 -x_2_5 0; x_1_4 0 -x_1_2 0 -x_3_4 0 x_4_5 x_2_3 0 -x_2_5; x_1_5 0 0 -x_1_2 -x_3_5 -x_4_5 0 0 x_2_3 x_2_4; 0 x_1_4 -x_1_3 0 x_2_4 -x_2_3 0 0 x_4_5 -x_3_5; 0 x_1_5 0 -x_1_3 x_2_5 0 -x_2_3 -x_4_5 0 x_3_4; 0 0 x_1_5 -x_1_4 0 x_2_5 -x_2_4 x_3_5 -x_3_4 0]"
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] == sparse(Int64[], Int64[], fmpq[], 1, 1)

                b = PseudographDeformBasis{fmpq}(sp, 0:3)
                @test length(collect(b)) == 4
                @test pbwdeforms_all(sp, b; special_return=SparseMatrixCSC)[1] ==
                      sparse([2, 3, 4], [2, 3, 4], fmpq[1, 1, 1], 4, 4)
            end
        end

    end

end
