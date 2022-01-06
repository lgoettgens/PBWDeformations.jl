C = DefaultScalarType

@testset ExtendedTestSet "All PBWDeformations.SmashProductLie tests" begin
    @testset "get_matrix_rep" begin
        @testset "A_$n" for n in 1:6
            matrixRep = PD.get_matrix_rep('A', n)
            @test length(matrixRep) == (n+1)^2-1
            @test all(mat -> size(mat) == (n+1, n+1), matrixRep)

            if n == 1
                @test matrixRep[1] in [c*[0 1; 0 0] for c in [-1, 1]]
                @test matrixRep[2] in [c*[0 0; 1 0] for c in [-1, 1]]
                @test matrixRep[3] in [c*[1 0; 0 -1] for c in [-1, 1]]
            elseif n == 2
                for i in 1:6
                    @test length([x for x in vcat(matrixRep[i]...) if x != 0]) == 1
                end
                @test matrixRep[7] in [c*[1 0 0; 0 -1 0; 0 0 0] for c in [-1, 1]]
                @test matrixRep[8] in [c*[0 0 0; 0 1 0; 0 0 -1] for c in [-1, 1]]
            end

        end

        @testset "B_$n" for n in 2:6
            matrixRep = PD.get_matrix_rep('B', n)
            @test length(matrixRep) == 2*n^2+n
            @test all(mat -> size(mat) == (2n+1, 2n+1), matrixRep)
        end

        @testset "C_$n" for n in 2:6
            matrixRep = PD.get_matrix_rep('C', n)
            @test length(matrixRep) == 2*n^2+n
            @test all(mat -> size(mat) == (2n, 2n), matrixRep)
        end

        @testset "D_$n" for n in 4:6
            matrixRep = PD.get_matrix_rep('D', n)
            @test length(matrixRep) == 2*n^2-n
            @test all(mat -> size(mat) == (2n, 2n), matrixRep)
        end

    end

end
