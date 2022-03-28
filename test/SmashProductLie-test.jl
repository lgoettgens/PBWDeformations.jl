@testset ExtendedTestSet "All SmashProductLie.jl tests" begin
    
    @testset "smash_product_lie constructor using GAP; R = $R" for R in [QQ, PolynomialRing(QQ, ["x","y","z"])[1], PolynomialRingSparse(QQ, ["x","y","z"])[1]]
        @testset "consistency for $(dynkin)_$n with hw $lambda" for (dynkin, n, lambda) in [('A', 1, [1]), ('A', 2, [1,1]), ('B', 2, [1,0])]
            sp, (baseL, baseV) = smash_product_lie(R, dynkin, n, lambda)
            
            @test sp.dimL == length(sp.baseL)
            @test sp.dimV == length(sp.baseV)
            @test sp.baseL == baseL
            @test sp.baseV == baseV
            @test sp.coeff_ring == R

            @test ngens(sp) == (sp.dimL, sp.dimV)
            @test gens(sp) == (sp.baseL, sp.baseV)

            showOutput = @test_nowarn sprint(show, sp)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            # @test_broken occursin(string(dynkin), showOutput)
            # @test_broken occursin(string(lambda), showOutput)
        end

        @testset "commutators for A_1 with hw [1]" begin
            sp, (baseL, baseV) = smash_product_lie(R, 'A', 1, [1])
            
            @test sp.dimL == 3
            @test sp.dimV == 2

            x = baseL[1]
            y = baseL[2]
            h = baseL[3]
            v1 = baseV[1]
            v2 = baseV[2]

            # sl_2 relations
            @test comm(x, x; strict=true) == 0
            @test comm(x, y; strict=true) == h
            @test comm(x, h; strict=true) == -2*x
            @test comm(y, x; strict=true) == -h
            @test comm(y, y; strict=true) == 0
            @test comm(y, h; strict=true) == 2*y
            @test comm(h, x; strict=true) == 2*x
            @test comm(h, y; strict=true) == -2*y
            @test comm(h, h; strict=true) == 0

            # natural representation relations
            @test comm(x, v1; strict=true) == 0
            @test comm(x, v2; strict=true) == v1
            @test comm(y, v1; strict=true) == v2
            @test comm(y, v2; strict=true) == 0
            @test comm(h, v1; strict=true) == v1
            @test comm(h, v2; strict=true) == -v2
        end

        @testset "commutators for A_2 with hw [1,1]" begin
            sp, (baseL, baseV) = smash_product_lie(R, 'A', 2, [1,0])
            
            @test sp.dimL == 8
            @test sp.dimV == 3

            x12 = baseL[1]
            x23 = -baseL[2] # minus due to representation in GAP 
            x13 = baseL[3]
            y12 = baseL[4]
            y23 = -baseL[5] # minus due to representation in GAP 
            y13 = baseL[6]
            h1 = baseL[7]
            h2 = baseL[8]
            v1 = baseV[1]
            v2 = baseV[2]
            v3 = baseV[3]

            # sl_3 relations
            @test comm(x12, x23; strict=true) == x13
            @test comm(x12, x13; strict=true) == 0
            @test comm(x23, x13; strict=true) == 0
            @test comm(y12, y23; strict=true) == -y13
            @test comm(y12, y13; strict=true) == 0
            @test comm(y23, y13; strict=true) == 0
            @test comm(h1, h2; strict=true) == 0
            @test comm(x12, y12; strict=true) == h1
            @test comm(x23, y23; strict=true) == h2
            @test comm(x13, y13; strict=true) == h1+h2
            @test comm(h1, x12; strict=true) == 2*x12
            @test comm(h1, x23; strict=true) == -x23
            @test comm(h1, x13; strict=true) == x13
            @test comm(h1, y12; strict=true) == -2*y12
            @test comm(h1, y23; strict=true) == y23
            @test comm(h1, y13; strict=true) == -y13
            @test comm(h2, x12; strict=true) == -x12
            @test comm(h2, x23; strict=true) == 2*x23
            @test comm(h2, x13; strict=true) == x13
            @test comm(h2, y12; strict=true) == y12
            @test comm(h2, y23; strict=true) == -2*y23
            @test comm(h2, y13; strict=true) == -y13


            # natural representation relations
            @test comm(x12, v1; strict=true) == 0
            @test comm(x23, v1; strict=true) == 0
            @test comm(x13, v1; strict=true) == 0
            @test comm(x12, v2; strict=true) == v1
            @test comm(x23, v2; strict=true) == 0
            @test comm(x13, v2; strict=true) == 0
            @test comm(x12, v3; strict=true) == 0
            @test comm(x23, v3; strict=true) == v2
            @test comm(x13, v3; strict=true) == v1
            @test comm(y12, v1; strict=true) == v2
            @test comm(y23, v1; strict=true) == 0
            @test comm(y13, v1; strict=true) == v3
            @test comm(y12, v2; strict=true) == 0
            @test comm(y23, v2; strict=true) == v3
            @test comm(y13, v2; strict=true) == 0
            @test comm(y12, v3; strict=true) == 0
            @test comm(y23, v3; strict=true) == 0
            @test comm(y13, v3; strict=true) == 0
            @test comm(h1, v1; strict=true) == v1
            @test comm(h2, v1; strict=true) == 0
            @test comm(h1, v2; strict=true) == -v2
            @test comm(h2, v2; strict=true) == v2
            @test comm(h1, v3; strict=true) == 0
            @test comm(h2, v3; strict=true) == -v3
        end

    end

    @testset "smash_product_lie constructor without GAP; R = $R" for R in [QQ, PolynomialRing(QQ, ["x","y","z"])[1], PolynomialRingSparse(QQ, ["x","y","z"])[1]]
        @testset "A_1 with hw [1]" begin
            dimL = 3
            dimV = 2
            symbL = ["a$i" for i in 1:dimL] # different to standard implementation
            symbV = ["b$i" for i in 1:dimV]
            struct_const_L = permutedims(reshape(Vector{Tuple{Int, Int}}[
                [],         [(1, 3)],   [(-2, 1)],
                [(-1, 3)],  [],         [(2, 2)],
                [(2, 1)],   [(-2, 2)],  [],
            ], 3, 3), (2, 1))
            struct_const_V = permutedims(reshape(Vector{Tuple{Int, Int}}[
                [],        [(1, 1)],
                [(1, 2)],  [],
                [(1, 1)],  [(-1, 2)],
            ], 2, 3), (2, 1))
            
            sp, (baseL, baseV) = smash_product_lie(R, symbL, symbV, struct_const_L, struct_const_V)
            
            @test dimL == sp.dimL == length(sp.baseL)
            @test dimV == sp.dimV == length(sp.baseV)
            @test sp.baseL == baseL
            @test sp.baseV == baseV
            @test sp.coeff_ring == R
            @test issetequal(sp.alg.S, map(Symbol, [symbL; symbV]))

            @test ngens(sp) == (sp.dimL, sp.dimV)
            @test gens(sp) == (sp.baseL, sp.baseV)

            showOutput = @test_nowarn sprint(show, sp)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))

            x = baseL[1]
            y = baseL[2]
            h = baseL[3]
            v1 = baseV[1]
            v2 = baseV[2]

            # sl_2 relations
            @test comm(x, x; strict=true) == 0
            @test comm(x, y; strict=true) == h
            @test comm(x, h; strict=true) == -2x
            @test comm(y, x; strict=true) == -h
            @test comm(y, y; strict=true) == 0
            @test comm(y, h; strict=true) == 2y
            @test comm(h, x; strict=true) == 2x
            @test comm(h, y; strict=true) == -2y
            @test comm(h, h; strict=true) == 0

            # natural representation relations
            @test comm(x, v1; strict=true) == 0
            @test comm(x, v2; strict=true) == v1
            @test comm(y, v1; strict=true) == v2
            @test comm(y, v2; strict=true) == 0
            @test comm(h, v1; strict=true) == v1
            @test comm(h, v2; strict=true) == -v2
        end


    end

       
end
