@testset ExtendedTestSet "All SmashProductLie.jl tests" begin
    @testset "smash_product; R = $R" for R in [QQ, PolynomialRing(QQ, ["x", "y", "z"])[1]]
        @testset "sl_2(QQ) ⋉ V" begin
            L = special_linear_lie_algebra(QQ, 2)
            V = standard_module(L)

            sp = smash_product(L, V)

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
            @test normal_form(comm(x, x), sp.rels) == 0
            @test normal_form(comm(x, y), sp.rels) == h
            @test normal_form(comm(x, h), sp.rels) == -2x
            @test normal_form(comm(y, x), sp.rels) == -h
            @test normal_form(comm(y, y), sp.rels) == 0
            @test normal_form(comm(y, h), sp.rels) == 2y
            @test normal_form(comm(h, x), sp.rels) == 2x
            @test normal_form(comm(h, y), sp.rels) == -2y
            @test normal_form(comm(h, h), sp.rels) == 0

            # natural representation relations
            @test normal_form(comm(x, v1), sp.rels) == 0
            @test normal_form(comm(x, v2), sp.rels) == v1
            @test normal_form(comm(y, v1), sp.rels) == v2
            @test normal_form(comm(y, v2), sp.rels) == 0
            @test normal_form(comm(h, v1), sp.rels) == v1
            @test normal_form(comm(h, v2), sp.rels) == -v2
        end

        @testset "so_4(QQ) ⋉ ⋀^2 V" begin
            L = special_orthogonal_lie_algebra(QQ, 4)
            V = exterior_power(standard_module(L), 2)

            sp = smash_product(L, V)

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
