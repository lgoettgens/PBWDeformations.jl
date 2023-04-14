@testset ExtendedTestSet "All LieAlgebraSymmetricPowerModule.jl tests" begin

    @testset "conformance tests for $desc" for (desc, L, V, parentT, elemT) in [
        (
            "S^3 V of sl_4(QQ)",
            special_linear_lie_algebra(QQ, 4),
            symmetric_power(standard_module(special_linear_lie_algebra(QQ, 4)), 3),
            LieAlgebraModule{fmpq},
            LieAlgebraModuleElem{fmpq},
        ),
        (
            "S^2 ⋀^2 V of so_4(QQ)",
            special_orthogonal_lie_algebra(QQ, 4),
            symmetric_power(exterior_power(standard_module(special_orthogonal_lie_algebra(QQ, 4)), 2), 2),
            LieAlgebraModule{fmpq},
            LieAlgebraModuleElem{fmpq},
        ),
        (
            "S^2 V of so_4(CL(4))",
            special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4),
            symmetric_power(standard_module(special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)), 3),
            LieAlgebraModule{nf_elem},
            LieAlgebraModuleElem{nf_elem},
        ),
    ]
        lie_algebra_module_conformance_test(L, V, parentT, elemT)
    end
end
