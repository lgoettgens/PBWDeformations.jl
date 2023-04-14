@testset ExtendedTestSet "All LieAlgebraAbstractModule.jl tests" begin

    R = QQ
    sc = Matrix{SRow{elem_type(R)}}(undef, 3, 2)
    sc[1, 1] = sparse_row(R, [1, 2], [0, 0])
    sc[1, 2] = sparse_row(R, [1, 2], [1, 0])
    sc[2, 1] = sparse_row(R, [1, 2], [0, 1])
    sc[2, 2] = sparse_row(R, [1, 2], [0, 0])
    sc[3, 1] = sparse_row(R, [1, 2], [1, 0])
    sc[3, 2] = sparse_row(R, [1, 2], [0, -1])

    @testset "conformance tests for $desc" for (desc, L, V, parentT, elemT) in [
        (
            "V of sl_2(QQ)",
            special_linear_lie_algebra(QQ, 2),
            abstract_module(special_linear_lie_algebra(QQ, 2), 2, sc),
            LieAlgebraModule{fmpq},
            LieAlgebraModuleElem{fmpq},
        ),
        (
            "λ = [1,1,0] of sl_4(QQ)",
            special_linear_lie_algebra(QQ, 4),
            highest_weight_module(special_linear_lie_algebra(QQ, 4), [1, 1, 0]),
            LieAlgebraModule{fmpq},
            LieAlgebraModuleElem{fmpq},
        ),
        # (
        #     "λ = [1,1,0,0] of so_4(QQ)",
        #     special_orthogonal_lie_algebra(QQ, 4),
        #     highest_weight_module(special_linear_lie_algebra(QQ, 4), [1, 1, 0, 0]),
        #     LieAlgebraModule{fmpq},
        #     LieAlgebraModuleElem{fmpq},
        # ),
        (
            "λ = [0,1] of A_2(QQ)",
            lie_algebra(QQ, ('A', 2)),
            highest_weight_module(lie_algebra(QQ, ('A', 2)), [0, 1]),
            LieAlgebraModule{fmpq},
            LieAlgebraModuleElem{fmpq},
        ),
        (
            "λ = [0,1,1] of B_3(QQ)",
            lie_algebra(QQ, ('B', 3)),
            highest_weight_module(lie_algebra(QQ, ('B', 3)), [0, 1, 1]),
            LieAlgebraModule{fmpq},
            LieAlgebraModuleElem{fmpq},
        ),
    ]

        lie_algebra_module_conformance_test(L, V, parentT, elemT)
    end
end
