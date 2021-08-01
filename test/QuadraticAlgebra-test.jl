@testset ExtendedTestSet "All PBWDeformations.QuadraticAlgebra tests" begin
    @testset "containment" begin
        basis = map(test, 1:10)
        relTable = Dict()
        alg = PD.QuadraticAlgebra{Rational{Int64}, Nothing}(basis, relTable, nothing) :: PD.QuadraticAlgebra{Rational{Int64}, Nothing}
        x = test

        @test basis[1] in alg
        @test prod(basis) in alg

        for _ in 1:numRandomTests
            @test randAlgebraElement(basis) in alg
        end

        @test !(x(0) in alg)
        @test !(x(11) in alg)
        @test !(x(-1) in alg)
        @test !(x(1)*(x(2) + x(34))*x(5,6) in alg)
    end


    @testset "normalForm for abstract cases" begin
        @testset "tensor algebra over V with dim V = $n" for n in dimRandomTests
            x = test
            basis = map(x, 1:n)
            relTable = Dict{Tuple{BasisElement{Rational{Int64}}, BasisElement{Rational{Int64}}}, AlgebraElement{Rational{Int64}}}()
            alg = PD.QuadraticAlgebra{Rational{Int64}, Nothing}(basis, relTable, nothing)

            showOutput = @test_nowarn sprint(show, alg)

            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:n)))
                ind2 = shuffle(rand(1:n, rand(1:n)))
                @test normalForm(alg, x([ind1; ind2])-x([ind2; ind1])) ≐ normalForm(alg, comm(x(ind1), x(ind2)))
            end
        end

        @testset "symmetric algebra over V with dim V = $n" for n in dimRandomTests
            x = test
            basis = map(x, 1:n)
            relTable = Dict([((x(i), x(j)), AlgebraElement{Rational{Int64}}(x(j,i))) for i in 1:n for j in 1:i-1])
            alg = PD.QuadraticAlgebra{Rational{Int64}, Nothing}(basis, relTable, nothing)

            showOutput = @test_nowarn sprint(show, alg)

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:n, rand(1:2n)))
                @test normalForm(alg, x(ind)) ≐ x(sort(ind))
            end

            for _ in 1:numRandomTests
                ind1 = shuffle(rand(1:n, rand(1:n)))
                ind2 = shuffle(rand(1:n, rand(1:n)))
                @test normalForm(alg, x([ind1; ind2])-x([ind2; ind1])) ≐ normalForm(alg, comm(x(ind1), x(ind2))) ≐ 0
            end
        end

        @testset "exterior algebra over V with dim V = $n" for n in dimRandomTests
            x = test
            basis = map(x, 1:n)
            relTable = Dict([
                [((x(i), x(j)), -x(j,i)) for i in 1:n for j in 1:i-1];
                [((x(i), x(i)), AlgebraElement{Rational{Int64}}(0)) for i in 1:n];
            ])
            alg = PD.QuadraticAlgebra{Rational{Int64}, Nothing}(basis, relTable, nothing)

            showOutput = @test_nowarn sprint(show, alg)

            for _ in 1:numRandomTests
                ind = shuffle(rand(1:n, rand(1:2n)))
                uniqInd = unique(ind)

                if length(unique(ind)) < length(ind)
                    @test normalForm(alg, x(ind)) ≐ 0
                end
                @test normalForm(alg, x(uniqInd)) ≐ levicivita(sortperm(uniqInd)) * x(sort(uniqInd))
            end

            for _ in 1:numRandomTests
                ind1 = unique(shuffle(rand(1:n, rand(1:n))))
                ind2 = unique(shuffle(rand(1:n, rand(1:n))))
                @test normalForm(alg, x([ind1; ind2])-x([ind2; ind1])) ≐ normalForm(alg, comm(x(ind1), x(ind2)))
            end
        end

    end

    @testset "normalForm for smashProductLie" begin
        @testset "A_2 with hw [1,1]" begin
            sp = PD.smashProductLie('A', 2, [1,1])

            showOutput = @test_nowarn sprint(show, sp)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            @test occursin("A", showOutput)
            @test occursin("[1,1]", showOutput) || occursin("[1, 1]", showOutput)

            # Lie elements commute as usual
            @test normalForm(sp, comm(lie(7), lie(1))) ≐ 2*lie(1)  # [h_1,x_1] = 2x_1
            @test normalForm(sp, comm(lie(7), lie(2))) ≐ -lie(2)   # [h_1,x_2] = -x_1
            @test normalForm(sp, comm(lie(7), lie(3))) ≐ lie(3)    # [h_1,x_3] =  x_1
            @test normalForm(sp, comm(lie(1), lie(4))) ≐ lie(7)    # [x_1,y_1] =  h_1

            # Module elements do not commute at all
            @test normalForm(sp, mod(1, 2)) ≐ mod(1, 2)
            @test normalForm(sp, mod(2, 1)) ≐ mod(2, 1)
            @test normalForm(sp, mod(5, 1)) ≐ mod(5, 1)
            @test normalForm(sp, mod(8, 4)) ≐ mod(8, 4)
            @test normalForm(sp, mod(1, 5, 8, 4, 2)) ≐ mod(1, 5, 8, 4, 2)

            # Application commutators
            @test normalForm(sp, comm(lie(1), mod(1))) ≐ 0
            @test normalForm(sp, comm(lie(2), mod(1))) ≐ 0
            @test normalForm(sp, comm(lie(3), mod(1))) ≐ 0
            @test normalForm(sp, comm(lie(4), mod(1))) ≐ mod(2)
            @test normalForm(sp, comm(lie(5), mod(1))) ≐ mod(3)
            @test normalForm(sp, comm(lie(6), mod(1))) ≐ mod(5)
            @test normalForm(sp, comm(lie(7), mod(1))) ≐ mod(1)
            @test normalForm(sp, comm(lie(8), mod(1))) ≐ mod(1)
            @test normalForm(sp, comm(lie(4), mod(3))) ≐ mod(4)
            @test normalForm(sp, comm(lie(5), mod(2))) ≐ mod(4) - mod(5)

            # Some more complicated
            @test normalForm(sp, lie(7, 1, 2, 3)) ≐ 2*lie(1, 2, 3) + lie(1, 2, 3, 7)
        end

    end

    @testset "normalForm for smashProductSymmDeformLie" begin
        @testset "A_2 with hw [1,1]" begin
            sp = PD.smashProductSymmDeformLie('A', 2, [1,1])

            showOutput = @test_nowarn sprint(show, sp)
            @test occursin("smash product", lowercase(showOutput))
            @test occursin("lie algebra", lowercase(showOutput))
            @test occursin("A", showOutput)
            @test occursin("[1,1]", showOutput) || occursin("[1, 1]", showOutput)
            @test occursin("symmetric deformation", lowercase(showOutput))

            # Lie elements commute as usual
            @test normalForm(sp, comm(lie(7), lie(1))) ≐ 2lie(1)   # [h_1,x_1] = 2x_1
            @test normalForm(sp, comm(lie(7), lie(2))) ≐ -lie(2)   # [h_1,x_2] = -x_1
            @test normalForm(sp, comm(lie(7), lie(3))) ≐ lie(3)    # [h_1,x_3] =  x_1
            @test normalForm(sp, comm(lie(1), lie(4))) ≐ lie(7)    # [x_1,y_1] =  h_1

            # Module elements do commute
            @test normalForm(sp, mod(1, 2)) ≐ mod(1, 2)
            @test normalForm(sp, mod(2, 1)) ≐ mod(1, 2)
            @test normalForm(sp, mod(5, 1)) ≐ mod(1, 5)
            @test normalForm(sp, mod(8, 4)) ≐ mod(4, 8)
            @test normalForm(sp, mod(1, 5, 8, 4, 2)) ≐ mod(1, 2, 4, 5, 8)
            

            # Application commutators
            @test normalForm(sp, comm(lie(1), mod(1))) ≐ 0
            @test normalForm(sp, comm(lie(2), mod(1))) ≐ 0
            @test normalForm(sp, comm(lie(3), mod(1))) ≐ 0
            @test normalForm(sp, comm(lie(4), mod(1))) ≐ mod(2)
            @test normalForm(sp, comm(lie(5), mod(1))) ≐ mod(3)
            @test normalForm(sp, comm(lie(6), mod(1))) ≐ mod(5)
            @test normalForm(sp, comm(lie(7), mod(1))) ≐ mod(1)
            @test normalForm(sp, comm(lie(8), mod(1))) ≐ mod(1)
            @test normalForm(sp, comm(lie(4), mod(3))) ≐ mod(4)
            @test normalForm(sp, comm(lie(5), mod(2))) ≐ mod(4) - mod(5)

            # Some more complicated
            @test normalForm(sp, lie(7, 1, 2, 3)) ≐ 2*lie(1, 2, 3) + lie(1, 2, 3, 7)
        end

    end

    @testset "normalForm for groupAlgebra" begin
        @testset "symmetric group S3" begin
            ga = PD.groupAlgebraSymmetricGroup(3)

            function perm(p::Vararg{Int64}) :: BasisElement
                str = filter(c -> !isspace(c), length(p) == 1 ? string("(", p[1], ")") : string(p))
                return grp(findfirst(isequal(str), ga.extraData.permRep))
            end

            @test normalForm(ga, perm(1,2)^2) ≐ perm()
            @test normalForm(ga, perm(1,3)^2) ≐ perm()
            @test normalForm(ga, perm(2,3)^2) ≐ perm()
            @test normalForm(ga, perm(1,2,3)^3) ≐ perm()
            @test normalForm(ga, perm(1,3,2)^3) ≐ perm()

            for _ in 1:numRandomTests
                ind = rand(1:factorial(3), rand(1:6))
                @test normalForm(ga, perm() * grp(ind)) ≐ normalForm(ga, grp(ind) * perm()) ≐ normalForm(ga, grp(ind))
            end

            @test normalForm(ga, perm(1,2) * perm(2,3)) ≐ perm(1,3,2)
            @test normalForm(ga, perm(2,3) * perm(1,2)) ≐ perm(1,2,3)
        end

        @testset "dicyclic group Q8" begin
            ga = PD.groupAlgebraDicyclicGroup(8)

            e = grp(findfirst(isequal("()"), ga.extraData.permRep))
            temp = filter(i -> grp(i) != e && normalForm(ga, grp(i)^2) ≐ e, 1:8)
            @test length(temp) == 1
            minus_e = grp(first(temp))
            for i in 1:8
                if !(grp(i) in [e, minus_e])
                    @test normalForm(ga, grp(i)^2) ≐ minus_e
                end
            end
        end
    end

end
