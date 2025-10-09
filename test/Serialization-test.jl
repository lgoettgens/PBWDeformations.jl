@testset verbose=true "Serialization.jl tests" begin
    test_save_load_roundtrip = PBWDeformations.test_save_load_roundtrip

    QQi = cyclotomic_field(4)[1]
    QQx = polynomial_ring(QQ, 5)[1]
    mktempdir() do path
        combinatorial_cases = [
            (ArcDiagramUndirected, arc_diagram(Undirected, "DBAB", "AD"), arc_diagram(Undirected, "ABBD", "AFFHHDKLLNNK")),
            (ArcDiagramDirected, arc_diagram(Directed, "dbAB", "aD"), arc_diagram(Directed, "aBbD", "AfFhHdKlLnNk")),
            (GlnGraph, GlnGraph(2, 2, [true, false, true, false], [(1, 4), (3, 2)]), GlnGraph(4, 4, [true, false, true, false, true, false, true, false], [(1, 8), (3, 6), (5, 2), (7, 4)])),
        ]

        @testset verbose=true "$type" for (type, x, y) in combinatorial_cases
            test_save_load_roundtrip(path, x) do loaded
                @test loaded == x
            end

            test_save_load_roundtrip(path, y) do loaded
                @test loaded == y
            end

            test_save_load_roundtrip(path, (x, y)) do loaded
                @test length(loaded) == 2
                @test loaded == (x, y)
            end

            test_save_load_roundtrip(path, [x, y]) do loaded
                @test length(loaded) == 2
                @test loaded == [x, y]
            end

            test_save_load_roundtrip(path, Dict("foo" => x, "bar" => y)) do loaded
                @test length(loaded) == 2
                @test issetequal(values(loaded), [x, y])
            end
        end

        @testset verbose=true "SmashProductLie" begin
            @testset "R = $R, LieR = $LieR" for (LieR, R) in [(QQ, QQ), (QQi, QQi), (QQ, QQi), (QQ, QQx)]
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
                        L = general_linear_lie_algebra(LieR, 2)
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
                        @test loaded == a
                        @assert loaded.simplified == a.simplified == true
                        @test data(loaded) == data(a)
                    end
                    simplify(a) # ensure that `a` is simplified
                    test_save_load_roundtrip(path, a) do loaded
                        @test parent(loaded) === parent(a)
                        @test loaded.simplified == a.simplified == true
                        @test data(loaded) == data(a)
                        @test loaded == a
                    end

                    test_save_load_roundtrip(path, gens(sp)) do loaded
                        @test length(loaded) == ngens(sp)
                        @test all(loaded[i] == gen(sp, i) for i in 1:ngens(sp))
                    end
                end
            end
        end

        @testset verbose=true "SmashProductLieDeform" begin
            @testset "R = $R" for R in [QQ, QQi]
                instances = [
                    begin
                        L = special_linear_lie_algebra(R, 2)
                        V = standard_module(L)
                        sp = smash_product(R, L, V)
                        d = symmetric_deformation(sp)
                    end,
                    begin
                        L = special_orthogonal_lie_algebra(R, 4, identity_matrix(R, 4))
                        V = exterior_power_obj(standard_module(L), 2)
                        sp = smash_product(R, L, V)

                        kappa = zero_matrix(sp, dim(V), dim(V))
                        kappa[1, 2] = gen(sp, 1, :L)
                        kappa[2, 1] = -kappa[1, 2]
                        kappa[3, 4] = gen(sp, 2, :L)
                        kappa[4, 3] = -kappa[3, 4]
                        d = deform(sp, kappa)
                    end,
                    begin
                        n = 3
                        L = general_linear_lie_algebra(R, n)
                        V = direct_sum(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                        kappa = zero_matrix(sp, dim(V), dim(V))
                        for i in 1:n
                            kappa[i, n+i] = 1
                            kappa[n+i, i] = -1
                        end
                        d = deform(sp, kappa)
                    end,
                    begin
                        L = general_linear_lie_algebra(R, 2)
                        V = tensor_product(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                        b = GlnGraphDeformBasis(sp, 4:4)
                        kappa = only(Iterators.filter(m -> occursin("((1, 1), (GlnGraph(2, 2, Bool[1, 0, 1, 0], [(1, 4), (3, 2)]), [2, 0], [2]))", string(lookup_params(m, b))), b))
                        d = deform(sp, kappa)
                    end,
                ]

                @testset "instance $i" for (i, d) in enumerate(instances)
                    test_save_load_roundtrip(path, d) do loaded
                        # nothing, cause `d === loaded` anyway
                    end

                    a = sum(i * gen(d, i) for i in 1:ngens(d)) + prod(gens(d)) + 5*prod(gen(d, i) for i in ngens(d):-2:1)
                    test_save_load_roundtrip(path, a) do loaded
                        @test parent(loaded) === parent(a)
                        @test loaded.simplified == a.simplified
                        @test data(loaded) == data(a)
                        @test loaded == a
                        @assert loaded.simplified == a.simplified == true
                        @test data(loaded) == data(a)
                    end
                    simplify(a) # ensure that `a` is simplified
                    test_save_load_roundtrip(path, a) do loaded
                        @test parent(loaded) === parent(a)
                        @test loaded.simplified == a.simplified == true
                        @test data(loaded) == data(a)
                        @test loaded == a
                    end

                    test_save_load_roundtrip(path, gens(d)) do loaded
                        @test length(loaded) == ngens(d)
                        @test all(loaded[i] == gen(d, i) for i in 1:ngens(d))
                    end
                end
            end
        end

        @testset verbose=true "StdDeformBasis" begin
            @testset "R = $R" for R in [QQ, QQi]
                instances = [
                    begin
                        L = special_linear_lie_algebra(R, 2)
                        V = standard_module(L)
                        sp = smash_product(R, L, V)
                        b = StdDeformBasis(sp, 0:3)
                    end,
                    begin
                        L = special_orthogonal_lie_algebra(R, 4, identity_matrix(R, 4))
                        V = exterior_power_obj(standard_module(L), 2)
                        sp = smash_product(R, L, V)
                        b = StdDeformBasis(sp, 0:2)
                    end,
                    begin
                        L = general_linear_lie_algebra(R, 3)
                        V = direct_sum(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                        b = StdDeformBasis(sp, 1:2)
                    end,
                    begin
                        L = general_linear_lie_algebra(R, 2)
                        V = tensor_product(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                        b = StdDeformBasis(sp, 1:1)
                    end,
                ]

                @testset "instance $i" for (i, b) in enumerate(instances)
                    test_save_load_roundtrip(path, b) do loaded
                        @test loaded.sp === b.sp
                        @test loaded.degs == b.degs
                        @test collect(loaded) == collect(b)
                    end
                end
            end
        end

        @testset verbose=true "ArcDiagDeformBasis" begin
            @testset "R = $R" for R in [QQ, QQi]
                instances = [
                    begin
                        L = special_orthogonal_lie_algebra(R, 4, identity_matrix(R, 4))
                        V = exterior_power_obj(standard_module(L), 2)
                        sp = smash_product(R, L, V)
                        b = ArcDiagDeformBasis(sp, 0:2)
                    end,
                    begin
                        L = general_linear_lie_algebra(R, 3)
                        V = direct_sum(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                        b = ArcDiagDeformBasis(sp, 1:2)
                    end,
                    begin
                        L = general_linear_lie_algebra(R, 2)
                        V = tensor_product(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                        b = ArcDiagDeformBasis(sp, 1:1)
                    end,
                ]

                @testset "instance $i" for (i, b) in enumerate(instances)
                    test_save_load_roundtrip(path, b) do loaded
                        @test loaded.sp === b.sp
                        @test loaded.degs == b.degs
                        @test collect(loaded) == collect(b)
                    end
                end
            end
        end

        @testset verbose=true "GlnGraphDeformBasis" begin
            @testset "R = $R" for R in [QQ, QQi]
                instances = [
                    begin
                        L = general_linear_lie_algebra(R, 3)
                        V = direct_sum(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                        b = GlnGraphDeformBasis(sp, 1:3)
                    end,
                    begin
                        L = general_linear_lie_algebra(R, 2)
                        V = tensor_product(standard_module(L), dual(standard_module(L)))
                        sp = smash_product(R, L, V)
                        b = GlnGraphDeformBasis(sp, 0:2)
                    end,
                ]

                @testset "instance $i" for (i, b) in enumerate(instances)
                    test_save_load_roundtrip(path, b) do loaded
                        @test loaded.sp === b.sp
                        @test loaded.degs == b.degs
                        @test collect(loaded) == collect(b)
                    end
                end
            end
        end
    end
end
