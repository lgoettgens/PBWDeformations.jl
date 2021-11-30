@testset ExtendedTestSet "All FreeAlgebra.jl tests" begin
    @testset "deepcopy_internal" begin
        F = PD.free_algebra(QQ, ["a"])
        A, = gens(F)

        B = deepcopy(A)

        @test parent(A) === parent(B)
        @test A.length == B.length
        @test A == B && A !== B
        @test A.monoms == B.monoms && A.monoms !== B.monoms
        @test A.coeffs == B.coeffs && A.coeffs !== B.coeffs
    end

    @testset "change_base_ring" begin
        F_Q = PD.free_algebra(QQ, ["a"])
        F_Z = PD.change_base_ring(ZZ, F_Q)

        @test base_ring(F_Q) != base_ring(F_Z)
        @test F_Q.S === F_Z.S
    end

    @testset "ismonomial" begin
        F = PD.free_algebra(QQ, ["a"])
        A, = gens(F)

        @test ismonomial(A)
        @test ismonomial(A^2)

        @test !ismonomial(2A)
        @test !ismonomial(A + 1)
    end

    @testset "iszero" begin
        F = PD.free_algebra(QQ, ["a"])
        A, = gens(F)

        @test iszero(0A)
        @test !iszero(A)
    end

    @testset "isone" begin
        F = PD.free_algebra(QQ, ["a"])
        A, = gens(F)

        @test isone(0A + 1)
        @test !isone(A)
    end

    @testset "isgen" begin
        F = PD.free_algebra(QQ, ["a"])
        A, = gens(F)

        @test isgen(A)
        
        @test !isgen(2A)
        @test !isgen(A^2)
        @test !isgen(A + 1)
    end

    @testset "Base.:-" begin
        F = PD.free_algebra(QQ, ["a", "b"])
        A, B = gens(F)

        @test -A !== -A
        @test A - A == zero(F)
        A - 1
    end
end