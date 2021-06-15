numRandomTests = 10
dimRandomTests = [3, 10, 25, 100]

@testset ExtendedTestSet "All PBWDeformations.AlgebraWithCommutators tests" begin
    @testset "normalForm for abstract cases" begin
        @testset "tensor algebra over V with dim V = $n" for n in dimRandomTests
            basis = [(:basis, i) for i in 1:n]
            commTable = Dict()
            alg = PD.AlgebraWithCommutators{Nothing}(basis, commTable, nothing)

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:n, rand(1:2n)))
                if length(ind) > 0
                    @test PD.normalForm(alg, alg.x(ind...)) == alg.x(ind...)
                end
            end
        end

        @testset "symmetric algebra over V with dim V = $n" for n in dimRandomTests
            basis = [(:symm, i) for i in 1:n]
            commTable = Dict([((basis[i], basis[j]), []) for i in 1:n for j in 1:i-1])
            alg = PD.AlgebraWithCommutators{Nothing}(basis, commTable, nothing)

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:n, rand(1:2n)))
                if length(ind) > 0
                    @test PD.normalForm(alg, alg.x(ind...)) == alg.x(sort(ind)...)
                end
            end
        end

    end

    @testset "normalForm for smashProductLie" begin
        @testset "A_2 with hw [1,1]" begin
            sp = PD.smashProductLie('A', 2, [1,1])

            # Lie elements commute as usual
            @test PD.normalForm(sp, sp.x(8+7, 8+1) - sp.x(8+1, 8+7)) == 2*sp.x(8+1)     # [h_1,x_1] = 2x_1
            @test PD.normalForm(sp, sp.x(8+7, 8+2) - sp.x(8+2, 8+7)) == (-1)*sp.x(8+2)  # [h_1,x_2] = -x_1
            @test PD.normalForm(sp, sp.x(8+7, 8+3) - sp.x(8+3, 8+7)) == sp.x(8+3)     # [h_1,x_3] =  x_1
            @test PD.normalForm(sp, sp.x(8+1, 8+4) - sp.x(8+4, 8+1)) == sp.x(8+7)     # [x_1,y_1] = h_1

            # Module elements do not commute at all
            @test PD.normalForm(sp, sp.x(1, 2)) == sp.x(1, 2)
            @test PD.normalForm(sp, sp.x(2, 1)) == sp.x(2, 1)
            @test PD.normalForm(sp, sp.x(5, 1)) == sp.x(5, 1)
            @test PD.normalForm(sp, sp.x(8, 4)) == sp.x(8, 4)
            @test PD.normalForm(sp, sp.x(1, 5, 8, 4, 2)) == sp.x(1, 5, 8, 4, 2)

            # Application commutators
            @test PD.normalForm(sp, sp.x(8+1, 1) - sp.x(1, 8+1)) == sympify(0)
            @test PD.normalForm(sp, sp.x(8+2, 1) - sp.x(1, 8+2)) == sympify(0)
            @test PD.normalForm(sp, sp.x(8+3, 1) - sp.x(1, 8+3)) == sympify(0)
            @test PD.normalForm(sp, sp.x(8+4, 1) - sp.x(1, 8+4)) == sp.x(2)
            @test PD.normalForm(sp, sp.x(8+5, 1) - sp.x(1, 8+5)) == sp.x(3)
            @test PD.normalForm(sp, sp.x(8+6, 1) - sp.x(1, 8+6)) == sp.x(5)
            @test PD.normalForm(sp, sp.x(8+7, 1) - sp.x(1, 8+7)) == sp.x(1)
            @test PD.normalForm(sp, sp.x(8+8, 1) - sp.x(1, 8+8)) == sp.x(1)
            @test PD.normalForm(sp, sp.x(8+4, 3) - sp.x(3, 8+4)) == sp.x(4)
            @test PD.normalForm(sp, sp.x(8+5, 2) - sp.x(2, 8+5)) == sp.x(4) - sp.x(5)

            # Some more complicated
            @test PD.normalForm(sp, sp.x(8+7, 8+1, 8+2, 8+3)) == 2*sp.x(8+1, 8+2, 8+3) + sp.x(8+1, 8+2, 8+3, 8+7)
        end

    end

    @testset "normalForm for smashProductSymmDeformLie" begin
        @testset "A_2 with hw [1,1]" begin
            sp = PD.smashProductSymmDeformLie('A', 2, [1,1])

            # Lie elements commute as usual
            @test PD.normalForm(sp, sp.x(8+7, 8+1) - sp.x(8+1, 8+7)) == 2*sp.x(8+1)     # [h_1,x_1] = 2x_1
            @test PD.normalForm(sp, sp.x(8+7, 8+2) - sp.x(8+2, 8+7)) == (-1)*sp.x(8+2)  # [h_1,x_2] = -x_1
            @test PD.normalForm(sp, sp.x(8+7, 8+3) - sp.x(8+3, 8+7)) == sp.x(8+3)     # [h_1,x_3] =  x_1
            @test PD.normalForm(sp, sp.x(8+1, 8+4) - sp.x(8+4, 8+1)) == sp.x(8+7)     # [x_1,y_1] = h_1

            # Module elements do commute
            @test PD.normalForm(sp, sp.x(1, 2)) == sp.x(1, 2)
            @test PD.normalForm(sp, sp.x(2, 1)) == sp.x(1, 2)
            @test PD.normalForm(sp, sp.x(5, 1)) == sp.x(1, 5)
            @test PD.normalForm(sp, sp.x(8, 4)) == sp.x(4, 8)
            @test PD.normalForm(sp, sp.x(1, 5, 8, 4, 2)) == sp.x(1, 2, 4, 5, 8)
            

            # Application commutators
            @test PD.normalForm(sp, sp.x(8+1, 1) - sp.x(1, 8+1)) == sympify(0)
            @test PD.normalForm(sp, sp.x(8+2, 1) - sp.x(1, 8+2)) == sympify(0)
            @test PD.normalForm(sp, sp.x(8+3, 1) - sp.x(1, 8+3)) == sympify(0)
            @test PD.normalForm(sp, sp.x(8+4, 1) - sp.x(1, 8+4)) == sp.x(2)
            @test PD.normalForm(sp, sp.x(8+5, 1) - sp.x(1, 8+5)) == sp.x(3)
            @test PD.normalForm(sp, sp.x(8+6, 1) - sp.x(1, 8+6)) == sp.x(5)
            @test PD.normalForm(sp, sp.x(8+7, 1) - sp.x(1, 8+7)) == sp.x(1)
            @test PD.normalForm(sp, sp.x(8+8, 1) - sp.x(1, 8+8)) == sp.x(1)
            @test PD.normalForm(sp, sp.x(8+4, 3) - sp.x(3, 8+4)) == sp.x(4)
            @test PD.normalForm(sp, sp.x(8+5, 2) - sp.x(2, 8+5)) == sp.x(4) - sp.x(5)
            
            # Some more complicated
            @test PD.normalForm(sp, sp.x(8+7, 8+1, 8+2, 8+3)) == 2*sp.x(8+1, 8+2, 8+3) + sp.x(8+1, 8+2, 8+3, 8+7)
        end

    end
       
end
