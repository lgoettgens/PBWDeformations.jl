@testset ExtendedTestSet "All SmashProductLie.jl tests" begin
    
    @testset "smash_product_lie constructor using GAP; R = $R" for R in [QQ, PolynomialRing(QQ, ["x","y","z"])[1]]
        @testset "$(dynkin)_$n with hw $lambda" for (dynkin, n, lambda) in [('A', 1, [1]), ('A', 2, [1,1]), ('B', 2, [1,0])]
            sp, (baseL, baseV) = PD.smash_product_lie(R, dynkin, n, lambda)
            
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
            sp, (baseL, baseV) = PD.smash_product_lie(R, 'A', 1, [1])
            
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

    @testset "smash_product_lie constructor without GAP; R = $R" for R in [QQ, PolynomialRing(QQ, ["x","y","z"])[1]]
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
            
            sp, (baseL, baseV) = PD.smash_product_lie(R, symbL, symbV, struct_const_L, struct_const_V)
            
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
