numRandomTests = 10
dimRandomTests = [3, 10, 25, 100]

@testset ExtendedTestSet "All PBWDeformations.QuadraticAlgebra tests" begin
    @testset "products and exponentiation" begin
        basis = [(:basis, i) for i in 1:10]
        relTable = Dict()
        alg = PD.QuadraticAlgebra{Nothing}(basis, relTable, nothing)
        x = alg.x

        @test PD.normalForm(alg, x(1,2)*x(3,4)) == x(1,2,3,4)
        @test PD.normalForm(alg, (x(1) + x(1,2))*x(3,4)) == x(1,3,4) + x(1,2,3,4)
        @test PD.normalForm(alg, x(1)*(x(2) + x(3,4))*x(5,6)) == x(1,2,5,6) + x(1,3,4,5,6)

        @test PD.normalForm(alg, x(1)^2*x(2)*x(3)^2) == x(1,1,2,3,3)
        @test PD.normalForm(alg, (x(1)+x(2))^2) == x(1,1) + x(1,2) + x(2,1) + x(2,2)
        @test PD.normalForm(alg, x(1)^5) == x(1,1,1,1,1)
    end

    @testset "normalForm for abstract cases" begin
        @testset "tensor algebra over V with dim V = $n" for n in dimRandomTests
            basis = [(:basis, i) for i in 1:n]
            relTable = Dict()
            alg = PD.QuadraticAlgebra{Nothing}(basis, relTable, nothing)
            x = alg.x

            showOutput = @test_nowarn sprint(show, alg)

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:n, rand(1:2n)))
                @test PD.normalForm(alg, prod(map(x, ind))) == PD.normalForm(alg, x(ind...)) == x(ind...)
            end

            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:n)))
                ind2 = shuffle(rand(1:n, rand(1:n)))
                @test PD.normalForm(alg, x(ind1..., ind2...)-x(ind2..., ind1...)) == PD.comm(alg, x(ind1...), x(ind2...))
            end
        end

        @testset "symmetric algebra over V with dim V = $n" for n in dimRandomTests
            basis = [(:symm, i) for i in 1:n]
            relTable = Dict([((basis[i], basis[j]), [(1, [basis[j], basis[i]])]) for i in 1:n for j in 1:i-1])
            alg = PD.QuadraticAlgebra{Nothing}(basis, relTable, nothing)
            x = alg.x

            showOutput = @test_nowarn sprint(show, alg)

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:n, rand(1:2n)))
                @test PD.normalForm(alg, x(ind...)) == x(sort(ind)...)
            end

            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:n)))
                ind2 = shuffle(rand(1:n, rand(1:n)))
                @test PD.normalForm(alg, x(ind1..., ind2...)-x(ind2..., ind1...)) == PD.comm(alg, x(ind1...), x(ind2...))
            end
        end

        @testset "exterior algebra over V with dim V = $n" for n in dimRandomTests
            basis = [(:alt, i) for i in 1:n]
            relTable = Dict([
                [((basis[i], basis[j]), [(-1, [basis[j], basis[i]])]) for i in 1:n for j in 1:i-1]...,
                [((basis[i], basis[i]), []) for i in 1:n]...,
            ])
            alg = PD.QuadraticAlgebra{Nothing}(basis, relTable, nothing)
            x = alg.x

            showOutput = @test_nowarn sprint(show, alg)

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:n, rand(1:2n)))
                uniqInd = unique(ind)

                if length(unique(ind)) < length(ind)
                    @test PD.normalForm(alg, x(ind...)) == sympify(0)
                end
                @test PD.normalForm(alg, x(uniqInd...)) == levicivita(sortperm(uniqInd)) * x(sort(uniqInd)...)
            end

            for _ in 1:numRandomTests
                ind1 = unique(shuffle(rand(1:n, rand(1:n))))
                ind2 = unique(shuffle(rand(1:n, rand(1:n))))
                @test PD.normalForm(alg, x(ind1..., ind2...)-x(ind2..., ind1...)) == PD.comm(alg, x(ind1...), x(ind2...))
            end
        end

    end

    @testset "normalForm for smashProductLie" begin
        @testset "A_2 with hw [1,1]" begin
            sp = PD.smashProductLie('A', 2, [1,1])
            x = sp.x

            showOutput = @test_nowarn sprint(show, sp)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            @test occursin("A", showOutput)
            @test occursin("[1,1]", showOutput) || occursin("[1, 1]", showOutput)

            # Lie elements commute as usual
            @test PD.normalForm(sp, x(8+7, 8+1) - x(8+1, 8+7)) == 2*x(8+1)     # [h_1,x_1] = 2x_1
            @test PD.normalForm(sp, x(8+7, 8+2) - x(8+2, 8+7)) == (-1)*x(8+2)  # [h_1,x_2] = -x_1
            @test PD.normalForm(sp, x(8+7, 8+3) - x(8+3, 8+7)) == x(8+3)     # [h_1,x_3] =  x_1
            @test PD.normalForm(sp, x(8+1, 8+4) - x(8+4, 8+1)) == x(8+7)     # [x_1,y_1] = h_1

            # Module elements do not commute at all
            @test PD.normalForm(sp, x(1, 2)) == x(1, 2)
            @test PD.normalForm(sp, x(2, 1)) == x(2, 1)
            @test PD.normalForm(sp, x(5, 1)) == x(5, 1)
            @test PD.normalForm(sp, x(8, 4)) == x(8, 4)
            @test PD.normalForm(sp, x(1, 5, 8, 4, 2)) == x(1, 5, 8, 4, 2)

            # Application commutators
            @test PD.normalForm(sp, x(8+1, 1) - x(1, 8+1)) == sympify(0)
            @test PD.normalForm(sp, x(8+2, 1) - x(1, 8+2)) == sympify(0)
            @test PD.normalForm(sp, x(8+3, 1) - x(1, 8+3)) == sympify(0)
            @test PD.normalForm(sp, x(8+4, 1) - x(1, 8+4)) == x(2)
            @test PD.normalForm(sp, x(8+5, 1) - x(1, 8+5)) == x(3)
            @test PD.normalForm(sp, x(8+6, 1) - x(1, 8+6)) == x(5)
            @test PD.normalForm(sp, x(8+7, 1) - x(1, 8+7)) == x(1)
            @test PD.normalForm(sp, x(8+8, 1) - x(1, 8+8)) == x(1)
            @test PD.normalForm(sp, x(8+4, 3) - x(3, 8+4)) == x(4)
            @test PD.normalForm(sp, x(8+5, 2) - x(2, 8+5)) == x(4) - x(5)

            # Some more complicated
            @test PD.normalForm(sp, x(8+7, 8+1, 8+2, 8+3)) == 2*x(8+1, 8+2, 8+3) + x(8+1, 8+2, 8+3, 8+7)

            n = length(sp.basis)
            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:6)))
                ind2 = shuffle(rand(1:n, rand(1:6)))
                @test PD.normalForm(sp, x(ind1..., ind2...)-x(ind2..., ind1...)) == PD.comm(sp, x(ind1...), x(ind2...))
            end
        end

    end

    @testset "normalForm for smashProductSymmDeformLie" begin
        @testset "A_2 with hw [1,1]" begin
            sp = PD.smashProductSymmDeformLie('A', 2, [1,1])
            x = sp.x

            showOutput = @test_nowarn sprint(show, sp)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            @test occursin("A", showOutput)
            @test occursin("[1,1]", showOutput) || occursin("[1, 1]", showOutput)
            @test occursin("symmetric deformation", lowercase(showOutput))

            # Lie elements commute as usual
            @test PD.normalForm(sp, x(8+7, 8+1) - x(8+1, 8+7)) == 2*x(8+1)     # [h_1,x_1] = 2x_1
            @test PD.normalForm(sp, x(8+7, 8+2) - x(8+2, 8+7)) == (-1)*x(8+2)  # [h_1,x_2] = -x_1
            @test PD.normalForm(sp, x(8+7, 8+3) - x(8+3, 8+7)) == x(8+3)     # [h_1,x_3] =  x_1
            @test PD.normalForm(sp, x(8+1, 8+4) - x(8+4, 8+1)) == x(8+7)     # [x_1,y_1] = h_1

            # Module elements do commute
            @test PD.normalForm(sp, x(1, 2)) == x(1, 2)
            @test PD.normalForm(sp, x(2, 1)) == x(1, 2)
            @test PD.normalForm(sp, x(5, 1)) == x(1, 5)
            @test PD.normalForm(sp, x(8, 4)) == x(4, 8)
            @test PD.normalForm(sp, x(1, 5, 8, 4, 2)) == x(1, 2, 4, 5, 8)
            

            # Application commutators
            @test PD.normalForm(sp, x(8+1, 1) - x(1, 8+1)) == sympify(0)
            @test PD.normalForm(sp, x(8+2, 1) - x(1, 8+2)) == sympify(0)
            @test PD.normalForm(sp, x(8+3, 1) - x(1, 8+3)) == sympify(0)
            @test PD.normalForm(sp, x(8+4, 1) - x(1, 8+4)) == x(2)
            @test PD.normalForm(sp, x(8+5, 1) - x(1, 8+5)) == x(3)
            @test PD.normalForm(sp, x(8+6, 1) - x(1, 8+6)) == x(5)
            @test PD.normalForm(sp, x(8+7, 1) - x(1, 8+7)) == x(1)
            @test PD.normalForm(sp, x(8+8, 1) - x(1, 8+8)) == x(1)
            @test PD.normalForm(sp, x(8+4, 3) - x(3, 8+4)) == x(4)
            @test PD.normalForm(sp, x(8+5, 2) - x(2, 8+5)) == x(4) - x(5)
            
            # Some more complicated
            @test PD.normalForm(sp, x(8+7, 8+1, 8+2, 8+3)) == 2*x(8+1, 8+2, 8+3) + x(8+1, 8+2, 8+3, 8+7)

            n = length(sp.basis)
            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:6)))
                ind2 = shuffle(rand(1:n, rand(1:6)))
                @test PD.normalForm(sp, x(ind1..., ind2...)-x(ind2..., ind1...)) == PD.comm(sp, x(ind1...), x(ind2...))
            end
        end

    end

end
