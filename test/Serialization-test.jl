@testset verbose=true "Serialization.jl tests" begin
    test_save_load_roundtrip = PBWDeformations.test_save_load_roundtrip

    QQi = cyclotomic_field(4)[1]
    mktempdir() do path
        @testset verbose=true "SmashProductLie" begin
            @testset "R = $R, LieR = $LieR" for (LieR, R) in [(QQ, QQ), (QQ, QQi), (QQi, QQi)]
                instances = [
                    begin
                        L = special_linear_lie_algebra(LieR, 2)
                        V = standard_module(L)
                        sp = smash_product(R, L, V)
                    end,
                    begin
                        L = special_orthogonal_lie_algebra(LieR, 4, identity_matrix(LieR, 4))
                        V = exterior_power_obj(standard_module(L), 2)
                        sp = smash_product(R, L, V)
                    end,
                    begin
                        L = general_linear_lie_algebra(LieR, 3)
                        V = direct_sum(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                    end,
                    begin
                        L = general_linear_lie_algebra(LieR, 3)
                        V = tensor_product(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                    end,
                ]

                @testset "instance $i" for (i, sp) in enumerate(instances)
                    test_save_load_roundtrip(path, sp) do loaded
                        # nothing, cause `sp === loaded` anyway
                    end

                    a = sum(i * gen(sp, i) for i in 1:ngens(sp)) + prod(gens(sp)) + 5*prod(gen(sp, i) for i in ngens(sp):-2:1)
                    test_save_load_roundtrip(path, a) do loaded
                        @test parent(loaded) === parent(a)
                        @test loaded.simplified == a.simplified
                        @test data(loaded) == data(a)
                        # don't call `loaded == a` as it would simplify both sides
                    end
                    simplify(a) # ensure that `a` is simplified
                    test_save_load_roundtrip(path, a) do loaded
                        @test parent(loaded) === parent(a)
                        @test loaded.simplified == a.simplified
                        @test data(loaded) == data(a)
                        # don't call `loaded == a` as it would simplify both sides
                    end

                    test_save_load_roundtrip(path, gens(sp)) do loaded
                        @test length(loaded) == ngens(sp)
                        @test all(loaded[i] == gen(sp, i) for i in 1:ngens(sp))
                    end
                end
            end
        end
    end
end
