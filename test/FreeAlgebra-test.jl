@testset ExtendedTestSet "All FreeAlgebra.jl tests" begin
    @testset "deepcopy_internal" begin
        A, (x,) = free_algebra(QQ, ["x"])

        y = deepcopy(x)

        @test parent(x) === parent(y)
        @test x.length == y.length
        @test x == y && x !== y
        @test x.monoms == y.monoms && x.monoms !== y.monoms
        @test x.coeffs == y.coeffs && x.coeffs !== y.coeffs
    end

    @testset "change_base_ring" begin
        A_Q, _ = free_algebra(QQ, ["x"])
        A_Z = change_base_ring(ZZ, A_Q)

        @test base_ring(A_Q) != base_ring(A_Z)
        @test A_Q.S === A_Z.S
    end

    @testset "is_monomial" begin
        A, (x,) = free_algebra(QQ, ["x"])

        @test is_monomial(x)
        @test is_monomial(x^2)

        @test !is_monomial(2x)
        @test !is_monomial(x + 1)
    end

    @testset "iszero" begin
        A, (x,) = free_algebra(QQ, ["x"])

        @test iszero(0 * x)
        @test iszero(x * 0)
        @test !iszero(x)
    end

    @testset "isone" begin
        A, (x,) = free_algebra(QQ, ["x"])

        @test isone(0 * x + 1)
        @test !isone(x)
    end

    @testset "is_gen" begin
        A, (x,) = free_algebra(QQ, ["x"])

        @test is_gen(x)

        @test !is_gen(2x)
        @test !is_gen(x^2)
        @test !is_gen(x + 1)
    end

    @testset "Base.:-" begin
        A, (x, y) = free_algebra(QQ, ["x", "y"])

        @test -x !== -x && -x == -x
        @test x - x == zero(A) == -x + x
        @test x + 1 - x == 1
        @test x - y == -y + x == x + (-1) * y == -1 * (y - x)
    end

    @testset "parent_type" begin
        A, (x,) = free_algebra(QQ, ["x"])

        # @test parent_type(A) == FreeAlgebra{fmpq}
    end

    @testset "show" begin
        A, (x,) = free_algebra(QQ, ["x"])

        showOutput = @test_nowarn sprint(show, A)
        @test occursin("free algebra", lowercase(showOutput))
    end
end
