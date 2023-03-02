@testset ExtendedTestSet "All SmashProductLie.jl tests" begin

    @testset "smash_product_lie_highest_weight constructor using GAP; R = $R" for R in [
        QQ,
        PolynomialRing(QQ, ["x", "y", "z"])[1],
    ]
        @testset "consistency for $(dynkin)_$n with hw $lambda" for (dynkin, n, lambda) in
                                                                    [('A', 1, [1]), ('A', 2, [1, 1]), ('B', 2, [1, 0])]
            sp, (basisL, basisV) = smash_product_lie_highest_weight(R, dynkin, n, lambda)

            @test sp.dimL == length(sp.basisL)
            @test sp.dimV == length(sp.basisV)
            @test sp.basisL == basisL
            @test sp.basisV == basisV
            @test sp.coeff_ring == R

            @test ngens(sp) == (sp.dimL, sp.dimV)
            @test gens(sp) == (sp.basisL, sp.basisV)

            @test sp.info.dynkin == dynkin
            @test sp.info.n == n
            @test sp.info.lambda == lambda
            @test sp.info.constructive_basis == false

            showOutput = @test_nowarn sprint(show, sp)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            # @test_broken occursin(string(dynkin), showOutput)
            # @test_broken occursin(string(lambda), showOutput)
        end

        @testset "commutators for A_1 with hw [1]" begin
            sp, (basisL, basisV) = smash_product_lie_highest_weight(R, 'A', 1, [1])

            @test sp.dimL == 3
            @test sp.dimV == 2

            x = basisL[1]
            y = basisL[2]
            h = basisL[3]
            v1 = basisV[1]
            v2 = basisV[2]

            # sl_2 relations
            @test normal_form(comm(x, x), sp.rels) == 0
            @test normal_form(comm(x, y), sp.rels) == h
            @test normal_form(comm(x, h), sp.rels) == -2 * x
            @test normal_form(comm(y, x), sp.rels) == -h
            @test normal_form(comm(y, y), sp.rels) == 0
            @test normal_form(comm(y, h), sp.rels) == 2 * y
            @test normal_form(comm(h, x), sp.rels) == 2 * x
            @test normal_form(comm(h, y), sp.rels) == -2 * y
            @test normal_form(comm(h, h), sp.rels) == 0

            # natural representation relations
            @test normal_form(comm(x, v1), sp.rels) == 0
            @test normal_form(comm(x, v2), sp.rels) == v1
            @test normal_form(comm(y, v1), sp.rels) == v2
            @test normal_form(comm(y, v2), sp.rels) == 0
            @test normal_form(comm(h, v1), sp.rels) == v1
            @test normal_form(comm(h, v2), sp.rels) == -v2
        end

        @testset "commutators for A_2 with hw [1,1]" begin
            sp, (basisL, basisV) = smash_product_lie_highest_weight(R, 'A', 2, [1, 0])

            @test sp.dimL == 8
            @test sp.dimV == 3

            x12 = basisL[1]
            x23 = -basisL[2] # minus due to representation in GAP
            x13 = basisL[3]
            y12 = basisL[4]
            y23 = -basisL[5] # minus due to representation in GAP
            y13 = basisL[6]
            h1 = basisL[7]
            h2 = basisL[8]
            v1 = basisV[1]
            v2 = basisV[2]
            v3 = basisV[3]

            # sl_3 relations
            @test normal_form(comm(x12, x23), sp.rels) == x13
            @test normal_form(comm(x12, x13), sp.rels) == 0
            @test normal_form(comm(x23, x13), sp.rels) == 0
            @test normal_form(comm(y12, y23), sp.rels) == -y13
            @test normal_form(comm(y12, y13), sp.rels) == 0
            @test normal_form(comm(y23, y13), sp.rels) == 0
            @test normal_form(comm(h1, h2), sp.rels) == 0
            @test normal_form(comm(x12, y12), sp.rels) == h1
            @test normal_form(comm(x23, y23), sp.rels) == h2
            @test normal_form(comm(x13, y13), sp.rels) == h1 + h2
            @test normal_form(comm(h1, x12), sp.rels) == 2 * x12
            @test normal_form(comm(h1, x23), sp.rels) == -x23
            @test normal_form(comm(h1, x13), sp.rels) == x13
            @test normal_form(comm(h1, y12), sp.rels) == -2 * y12
            @test normal_form(comm(h1, y23), sp.rels) == y23
            @test normal_form(comm(h1, y13), sp.rels) == -y13
            @test normal_form(comm(h2, x12), sp.rels) == -x12
            @test normal_form(comm(h2, x23), sp.rels) == 2 * x23
            @test normal_form(comm(h2, x13), sp.rels) == x13
            @test normal_form(comm(h2, y12), sp.rels) == y12
            @test normal_form(comm(h2, y23), sp.rels) == -2 * y23
            @test normal_form(comm(h2, y13), sp.rels) == -y13


            # natural representation relations
            @test normal_form(comm(x12, v1), sp.rels) == 0
            @test normal_form(comm(x23, v1), sp.rels) == 0
            @test normal_form(comm(x13, v1), sp.rels) == 0
            @test normal_form(comm(x12, v2), sp.rels) == v1
            @test normal_form(comm(x23, v2), sp.rels) == 0
            @test normal_form(comm(x13, v2), sp.rels) == 0
            @test normal_form(comm(x12, v3), sp.rels) == 0
            @test normal_form(comm(x23, v3), sp.rels) == v2
            @test normal_form(comm(x13, v3), sp.rels) == v1
            @test normal_form(comm(y12, v1), sp.rels) == v2
            @test normal_form(comm(y23, v1), sp.rels) == 0
            @test normal_form(comm(y13, v1), sp.rels) == v3
            @test normal_form(comm(y12, v2), sp.rels) == 0
            @test normal_form(comm(y23, v2), sp.rels) == v3
            @test normal_form(comm(y13, v2), sp.rels) == 0
            @test normal_form(comm(y12, v3), sp.rels) == 0
            @test normal_form(comm(y23, v3), sp.rels) == 0
            @test normal_form(comm(y13, v3), sp.rels) == 0
            @test normal_form(comm(h1, v1), sp.rels) == v1
            @test normal_form(comm(h2, v1), sp.rels) == 0
            @test normal_form(comm(h1, v2), sp.rels) == -v2
            @test normal_form(comm(h2, v2), sp.rels) == v2
            @test normal_form(comm(h1, v3), sp.rels) == 0
            @test normal_form(comm(h2, v3), sp.rels) == -v3
        end

    end

    @testset "smash_product_lie constructor without GAP; R = $R" for R in [QQ, PolynomialRing(QQ, ["x", "y", "z"])[1]]
        @testset "A_1 with hw [1]" begin
            dimL = 3
            dimV = 2
            symbL = ["a$i" for i in 1:dimL] # different to standard implementation
            symbV = ["b$i" for i in 1:dimV]
            struct_const_L = permutedims(
                reshape(
                    Vector{Tuple{elem_type(R), Int}}[
                        [],
                        [(R(1), 3)],
                        [(R(-2), 1)],
                        [(R(-1), 3)],
                        [],
                        [(R(2), 2)],
                        [(R(2), 1)],
                        [(R(-2), 2)],
                        [],
                    ],
                    3,
                    3,
                ),
                (2, 1),
            )
            struct_const_V = permutedims(
                reshape(
                    Vector{Tuple{elem_type(R), Int}}[[], [(R(1), 1)], [(R(1), 2)], [], [(R(1), 1)], [(R(-1), 2)]],
                    2,
                    3,
                ),
                (2, 1),
            )

            sp, (basisL, basisV) = smash_product_lie(R, symbL, symbV, struct_const_L, struct_const_V)

            @test dimL == sp.dimL == length(sp.basisL)
            @test dimV == sp.dimV == length(sp.basisV)
            @test sp.basisL == basisL
            @test sp.basisV == basisV
            @test sp.coeff_ring == R
            @test issetequal(sp.alg.S, map(Symbol, [symbL; symbV]))

            @test ngens(sp) == (sp.dimL, sp.dimV)
            @test gens(sp) == (sp.basisL, sp.basisV)

            showOutput = @test_nowarn sprint(show, sp)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))

            x = basisL[1]
            y = basisL[2]
            h = basisL[3]
            v1 = basisV[1]
            v2 = basisV[2]

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

    end

end
