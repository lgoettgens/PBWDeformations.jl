@testset ExtendedTestSet "All FreeAlgebra.jl tests" begin
    @testset "deepcopy_internal" begin
        F, (A,) = free_algebra(QQ, ["a"])

        B = deepcopy(A)

        @test parent(A) === parent(B)
        @test A.length == B.length
        @test A == B && A !== B
        @test A.monoms == B.monoms && A.monoms !== B.monoms
        @test A.coeffs == B.coeffs && A.coeffs !== B.coeffs
    end

    @testset "change_base_ring" begin
        F_Q, = free_algebra(QQ, ["a"])
        F_Z = change_base_ring(ZZ, F_Q)

        @test base_ring(F_Q) != base_ring(F_Z)
        @test F_Q.S === F_Z.S
    end

    @testset "is_monomial" begin
        F, (A,) = free_algebra(QQ, ["a"])

        @test is_monomial(A)
        @test is_monomial(A^2)

        @test !is_monomial(2A)
        @test !is_monomial(A + 1)
    end

    @testset "iszero" begin
        F, (A,) = free_algebra(QQ, ["a"])

        @test iszero(0A)
        @test !iszero(A)
    end

    @testset "isone" begin
        F, (A,) = free_algebra(QQ, ["a"])

        @test isone(0A + 1)
        @test !isone(A)
    end

    @testset "is_gen" begin
        F, (A,) = free_algebra(QQ, ["a"])

        @test is_gen(A)

        @test !is_gen(2A)
        @test !is_gen(A^2)
        @test !is_gen(A + 1)
    end

    @testset "Base.:-" begin
        F, (A, B) = free_algebra(QQ, ["a", "b"])

        @test -A !== -A && -A == -A
        @test A - A == zero(F) == -A + A
        @test A + 1 - A == 1
    end

    @testset "parent_type" begin
        F, (A,) = free_algebra(QQ, ["a"])

        # @test parent_type(A) == FreeAlgebra{fmpq}
    end
end
