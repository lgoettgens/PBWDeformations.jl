@testset ExtendedTestSet "All LieAlgebraStdModule.jl tests" begin

    @testset "constructors for R=$R, n=$n" for R in [QQ, cyclotomic_field(4)[1]], n in 1:5
        L = general_linear_lie_algebra(R, n)
        V = standard_module(L)
        @test dim(V) == n

        L = special_linear_lie_algebra(R, n)
        V = standard_module(L)
        @test dim(V) == n

        L = special_orthogonal_lie_algebra(R, n)
        V = standard_module(L)
        @test dim(V) == n
    end

    @testset "conformance tests for $desc" for (desc, L, V, parentT, elemT) in [
        (
            "V of sl_4(QQ)",
            special_linear_lie_algebra(QQ, 4),
            standard_module(special_linear_lie_algebra(QQ, 4)),
            LieAlgebraModule{fmpq},
            LieAlgebraModuleElem{fmpq},
        ),
        (
            "V of so_4(QQ)",
            special_orthogonal_lie_algebra(QQ, 4),
            standard_module(special_orthogonal_lie_algebra(QQ, 4)),
            LieAlgebraModule{fmpq},
            LieAlgebraModuleElem{fmpq},
        ),
        (
            "V of sl_4(CF(4))",
            special_linear_lie_algebra(cyclotomic_field(4)[1], 4),
            standard_module(special_linear_lie_algebra(cyclotomic_field(4)[1], 4)),
            LieAlgebraModule{nf_elem},
            LieAlgebraModuleElem{nf_elem},
        ),
        (
            "V of so_4(CF(4))",
            special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4),
            standard_module(special_orthogonal_lie_algebra(cyclotomic_field(4)[1], 4)),
            LieAlgebraModule{nf_elem},
            LieAlgebraModuleElem{nf_elem},
        ),
    ]
        lie_algebra_module_conformance_test(L, V, parentT, elemT)
    end
end
