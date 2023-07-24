@testset "SmashProductLie.jl tests" begin
    @testset "smash_product; R = $R" for R in [QQ, polynomial_ring(QQ, ["x", "y", "z"])[1]]
        @testset "sl_2(QQ) ⋉ V" begin
            L = special_linear_lie_algebra(QQ, 2)
            V = standard_module(L)

            sp = smash_product(R, L, V)

            @test (@inferred coefficient_ring(sp)) == R
            @test (@inferred base_lie_algebra(sp)) == L
            @test (@inferred base_module(sp)) == V

            @test dim(L) == ngens(sp, :L)
            @test dim(L) == length(gens(sp, :L))
            @test dim(L) == 3

            @test dim(V) == ngens(sp, :V)
            @test dim(V) == length(gens(sp, :V))
            @test dim(V) == 2

            @test dim(L) + dim(V) == ngens(sp)
            @test dim(L) + dim(V) == length(gens(sp))
            @test dim(L) + dim(V) == 5

            x = gen(sp, 1, :L)
            y = gen(sp, 2, :L)
            h = gen(sp, 3, :L)
            v1 = gen(sp, 1, :V)
            v2 = gen(sp, 2, :V)

            # sl_2 relations
            @test comm(x, x) == 0
            @test comm(x, y) == h
            @test comm(x, h) == -2x
            @test comm(y, x) == -h
            @test comm(y, y) == 0
            @test comm(y, h) == 2y
            @test comm(h, x) == 2x
            @test comm(h, y) == -2y
            @test comm(h, h) == 0

            # natural representation relations
            @test comm(x, v1) == 0
            @test comm(x, v2) == v1
            @test comm(y, v1) == v2
            @test comm(y, v2) == 0
            @test comm(h, v1) == v1
            @test comm(h, v2) == -v2
        end

        @testset "so_4(QQ) ⋉ ⋀^2 V" begin
            L = special_orthogonal_lie_algebra(QQ, 4)
            V = exterior_power(standard_module(L), 2)

            sp = smash_product(R, L, V)

            @test (@inferred coefficient_ring(sp)) == R
            @test (@inferred base_lie_algebra(sp)) == L
            @test (@inferred base_module(sp)) == V

            @test dim(L) == ngens(sp, :L)
            @test dim(L) == length(gens(sp, :L))
            @test dim(L) == 6

            @test dim(V) == ngens(sp, :V)
            @test dim(V) == length(gens(sp, :V))
            @test dim(V) == 6

            @test dim(L) + dim(V) == ngens(sp)
            @test dim(L) + dim(V) == length(gens(sp))
            @test dim(L) + dim(V) == 12
        end

    end

end
