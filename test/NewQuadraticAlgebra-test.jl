@testset ExtendedTestSet "All PBWDeformations.NewQuadraticAlgebra tests" begin
    @testset "containment" begin
        eval(PD.createBasisFunctions(:b))
        relTable = Dict()
        alg = PD.QuadraticAlgebra{Nothing}(basis(collect(1:10)), relTable, nothing) :: PD.QuadraticAlgebra{Nothing}

        @test basis(1) in alg
        @test basis(collect(1:10)) in alg # product over all basis elements
        @test !(basis(0) in alg)
        @test !(basis(-1) in alg)
        @test !(basis(11) in alg)
        @test !(basis(1)*(basis(2) + basis(34))*basis(5,6) in alg)

        for _ in 1:numRandomTests
            @test randAlgebraElement(basis(collect(1:10))) in alg
        end
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
                @test normalForm(alg, prod(map(x, ind))) == normalForm(alg, x(ind...)) == x(ind...)
            end

            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:n)))
                ind2 = shuffle(rand(1:n, rand(1:n)))
                @test normalForm(alg, x(ind1..., ind2...)-x(ind2..., ind1...)) == normalForm(alg, PD.comm(x(ind1...), x(ind2...)))
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
                @test normalForm(alg, x(ind...)) == x(sort(ind)...)
            end

            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:n)))
                ind2 = shuffle(rand(1:n, rand(1:n)))
                @test normalForm(alg, x(ind1..., ind2...)-x(ind2..., ind1...)) == normalForm(alg, PD.comm(x(ind1...), x(ind2...)))
            end
        end

        @testset "exterior algebra over V with dim V = $n" for n in dimRandomTests
            basis = [(:alt, i) for i in 1:n]
            relTable = Dict([
                [((basis[i], basis[j]), [(-1, [basis[j], basis[i]])]) for i in 1:n for j in 1:i-1];
                [((basis[i], basis[i]), []) for i in 1:n];
            ])
            alg = PD.QuadraticAlgebra{Nothing}(basis, relTable, nothing)
            x = alg.x

            showOutput = @test_nowarn sprint(show, alg)

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:n, rand(1:2n)))
                uniqInd = unique(ind)

                if length(unique(ind)) < length(ind)
                    @test normalForm(alg, x(ind...)) == sympify(0)
                end
                @test normalForm(alg, x(uniqInd...)) == levicivita(sortperm(uniqInd)) * x(sort(uniqInd)...)
            end

            for _ in 1:numRandomTests
                ind1 = unique(shuffle(rand(1:n, rand(1:n))))
                ind2 = unique(shuffle(rand(1:n, rand(1:n))))
                @test normalForm(alg, x(ind1..., ind2...)-x(ind2..., ind1...)) == normalForm(alg, PD.comm(x(ind1...), x(ind2...)))
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
            @test normalForm(sp, x(8+7, 8+1) - x(8+1, 8+7)) == 2*x(8+1)     # [h_1,x_1] = 2x_1
            @test normalForm(sp, x(8+7, 8+2) - x(8+2, 8+7)) == (-1)*x(8+2)  # [h_1,x_2] = -x_1
            @test normalForm(sp, x(8+7, 8+3) - x(8+3, 8+7)) == x(8+3)     # [h_1,x_3] =  x_1
            @test normalForm(sp, x(8+1, 8+4) - x(8+4, 8+1)) == x(8+7)     # [x_1,y_1] = h_1

            # Module elements do not commute at all
            @test normalForm(sp, x(1, 2)) == x(1, 2)
            @test normalForm(sp, x(2, 1)) == x(2, 1)
            @test normalForm(sp, x(5, 1)) == x(5, 1)
            @test normalForm(sp, x(8, 4)) == x(8, 4)
            @test normalForm(sp, x(1, 5, 8, 4, 2)) == x(1, 5, 8, 4, 2)

            # Application commutators
            @test normalForm(sp, x(8+1, 1) - x(1, 8+1)) == sympify(0)
            @test normalForm(sp, x(8+2, 1) - x(1, 8+2)) == sympify(0)
            @test normalForm(sp, x(8+3, 1) - x(1, 8+3)) == sympify(0)
            @test normalForm(sp, x(8+4, 1) - x(1, 8+4)) == x(2)
            @test normalForm(sp, x(8+5, 1) - x(1, 8+5)) == x(3)
            @test normalForm(sp, x(8+6, 1) - x(1, 8+6)) == x(5)
            @test normalForm(sp, x(8+7, 1) - x(1, 8+7)) == x(1)
            @test normalForm(sp, x(8+8, 1) - x(1, 8+8)) == x(1)
            @test normalForm(sp, x(8+4, 3) - x(3, 8+4)) == x(4)
            @test normalForm(sp, x(8+5, 2) - x(2, 8+5)) == x(4) - x(5)

            # Some more complicated
            @test normalForm(sp, x(8+7, 8+1, 8+2, 8+3)) == 2*x(8+1, 8+2, 8+3) + x(8+1, 8+2, 8+3, 8+7)

            n = length(sp.basis)
            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:6)))
                ind2 = shuffle(rand(1:n, rand(1:6)))
                @test normalForm(sp, x(ind1..., ind2...)-x(ind2..., ind1...)) == normalForm(sp, PD.comm(x(ind1...), x(ind2...)))
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
            @test normalForm(sp, x(8+7, 8+1) - x(8+1, 8+7)) == 2*x(8+1)     # [h_1,x_1] = 2x_1
            @test normalForm(sp, x(8+7, 8+2) - x(8+2, 8+7)) == (-1)*x(8+2)  # [h_1,x_2] = -x_1
            @test normalForm(sp, x(8+7, 8+3) - x(8+3, 8+7)) == x(8+3)     # [h_1,x_3] =  x_1
            @test normalForm(sp, x(8+1, 8+4) - x(8+4, 8+1)) == x(8+7)     # [x_1,y_1] = h_1

            # Module elements do commute
            @test normalForm(sp, x(1, 2)) == x(1, 2)
            @test normalForm(sp, x(2, 1)) == x(1, 2)
            @test normalForm(sp, x(5, 1)) == x(1, 5)
            @test normalForm(sp, x(8, 4)) == x(4, 8)
            @test normalForm(sp, x(1, 5, 8, 4, 2)) == x(1, 2, 4, 5, 8)
            

            # Application commutators
            @test normalForm(sp, x(8+1, 1) - x(1, 8+1)) == sympify(0)
            @test normalForm(sp, x(8+2, 1) - x(1, 8+2)) == sympify(0)
            @test normalForm(sp, x(8+3, 1) - x(1, 8+3)) == sympify(0)
            @test normalForm(sp, x(8+4, 1) - x(1, 8+4)) == x(2)
            @test normalForm(sp, x(8+5, 1) - x(1, 8+5)) == x(3)
            @test normalForm(sp, x(8+6, 1) - x(1, 8+6)) == x(5)
            @test normalForm(sp, x(8+7, 1) - x(1, 8+7)) == x(1)
            @test normalForm(sp, x(8+8, 1) - x(1, 8+8)) == x(1)
            @test normalForm(sp, x(8+4, 3) - x(3, 8+4)) == x(4)
            @test normalForm(sp, x(8+5, 2) - x(2, 8+5)) == x(4) - x(5)
            
            # Some more complicated
            @test normalForm(sp, x(8+7, 8+1, 8+2, 8+3)) == 2*x(8+1, 8+2, 8+3) + x(8+1, 8+2, 8+3, 8+7)

            n = length(sp.basis)
            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:6)))
                ind2 = shuffle(rand(1:n, rand(1:6)))
                @test normalForm(sp, x(ind1..., ind2...)-x(ind2..., ind1...)) == normalForm(sp, PD.comm(x(ind1...), x(ind2...)))
            end
        end

    end

    @testset "normalForm for groupAlgebra" begin
        @testset "symmetric group S3" begin
            ga = PD.groupAlgebraSymmetricGroup(3)

            x = ga.x
            fromPerm(p) = findfirst(isequal(p), ga.extraData.permRep)
            @test normalForm(ga, x(fromPerm("(1,2)"))*x(fromPerm("(1,2)"))) == normalForm(ga, x(fromPerm("()")))
            @test normalForm(ga, x(fromPerm("(1,3)"))*x(fromPerm("(1,3)"))) == normalForm(ga, x(fromPerm("()")))
            @test normalForm(ga, x(fromPerm("(2,3)"))*x(fromPerm("(2,3)"))) == normalForm(ga, x(fromPerm("()")))
            @test normalForm(ga, x(fromPerm("(1,2,3)"))*x(fromPerm("(1,2,3)"))*x(fromPerm("(1,2,3)"))) == normalForm(ga, x(fromPerm("()")))
            @test normalForm(ga, x(fromPerm("(1,3,2)"))*x(fromPerm("(1,3,2)"))*x(fromPerm("(1,3,2)"))) == normalForm(ga, x(fromPerm("()")))
            for _ in 1:numRandomTests
                ind = rand(1:factorial(3), rand(1:6))
                @test normalForm(ga, x(fromPerm("()"))*x(ind...)) == normalForm(ga, x(ind...) * x(fromPerm("()"))) == normalForm(ga, x(ind...))
            end

            @test normalForm(ga, x(fromPerm("(1,2)"))*x(fromPerm("(2,3)"))) == normalForm(ga, x(fromPerm("(1,3,2)")))
            @test normalForm(ga, x(fromPerm("(2,3)"))*x(fromPerm("(1,2)"))) == normalForm(ga, x(fromPerm("(1,2,3)")))
        end

        @testset "dicyclic group Q8" begin
            ga = PD.groupAlgebraDicyclicGroup(8)

            x = ga.x
            fromPerm(p) = findfirst(isequal(p), ga.extraData.permRep)

            e = fromPerm("()")
            temp = filter(g -> g != e && normalForm(ga, x(g)*x(g)) == x(e), 1:8)
            @test length(temp) == 1
            minus_e = first(temp)
            for g in 1:8
                if g != e && g != minus_e
                    @test normalForm(ga, x(g)*x(g)) == x(minus_e)
                end
            end
        end
    end

end
