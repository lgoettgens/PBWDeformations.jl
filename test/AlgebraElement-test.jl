x = test
randNum() = rand(-20:20)
randNums(quantity) = rand(-20:20, quantity)
randLength(start=0) = rand(start:10)
randBasisElement() = x(randNum()) :: BasisElement
basis = x(collect(-20:20))
randMonomial() = randMonomial(basis)
randAlgebraElement() = randAlgebraElement(basis)

@testset ExtendedTestSet "All AlgebraElement.jl tests" begin
    @testset "test conversions" begin
        @test algebraElement() == algebraElement(0) == algebraElement(0) == AlgebraElement([]) == AlgebraElement()
        @test algebraElement(algebraElement(0)) == algebraElement(0)

        @test algebraElement(1) == algebraElement(1) == algebraElement(Monomial{BasisElement}()) == [(Coefficient(1), Monomial{BasisElement}())] :: AlgebraElement
        @test algebraElement(algebraElement(1)) == algebraElement(1)

        @test algebraElement(test(1)) == algebraElement([x(1)]) == [(Coefficient(1), [x(1)])] :: AlgebraElement
        @test algebraElement(algebraElement(x(1))) == algebraElement(x(1))

        @test algebraElement(x(1,2,3)) == [(Coefficient(1), x(1,2,3))] :: AlgebraElement
        @test algebraElement(algebraElement(x(1,2,3))) == algebraElement(x(1,2,3))

        a = [(Coefficient(1//3), [x(1)]), (Coefficient(-1), x(1,1))] :: AlgebraElement
        @test algebraElement(a) == a
    end

    @testset "test misc operations" begin
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test PD.collectSummands(map(algebraElement, nums)) == algebraElement(sum(nums))
        end

        for _ in 1:numRandomTests
            b = randBasisElement()
            @test PD.collectSummands(algebraElement(b)) == algebraElement(b)
        end

        for _ in 1:numRandomTests
            mon = randMonomial()
            @test PD.collectSummands(algebraElement(mon)) == algebraElement(mon)
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test a ≐ PD.collectSummands(a)
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test a ≐ shuffle(a)
        end

        for _ in 1:numRandomTests
            c1, c2, c3 = map(Coefficient, randNums(3))
            mon = randMonomial()
            mon2 = randMonomial()
            @test [(c1, mon), (c2, mon), (c3, mon2)] ≐ [(c3, mon2), (c1+c2, mon)]
        end

        for _ in 1:numRandomTests
            c1, c2 = map(Coefficient, randNums(2))
            mon = randMonomial()
            mon2 = randMonomial()
            @test [(c1, mon), (c2, mon2), (-c1, mon)] ≐ [(c2, mon2)]
        end
    end

    @testset "test addition" begin
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test sum(map(algebraElement, nums)) ≐ algebraElement(sum(nums))
        end

        for _ in 1:numRandomTests
            as = [randAlgebraElement() for _ in 1:randLength(2)]
            @test sum(shuffle(as)) ≐ sum(as)
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test 0 + a ≐ a + 0 ≐ a
        end

        for _ in 1:numRandomTests
            c = randNum()
            a = randAlgebraElement()
            @test c + a ≐ Coefficient(c) + a
        end

        for _ in 1:numRandomTests
            c = randNum()
            b = randBasisElement()
            m = randMonomial()
            a = randAlgebraElement()
            @test c + b + a + m ≐ sum(algebraElement, shuffle([c, b, a, m]))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            b = randBasisElement()
            @test [(Coefficient(l), [b])] ≐ sum(fill(b, l))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            m = randMonomial()
            @test [(Coefficient(l), m)] ≐ sum(fill(m, l))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            a = randAlgebraElement()
            b = algebraElement([(Coefficient(l*coeff), mon) for (coeff, mon) in a])
            @test sameSum(sum(fill(a,l)), b)
        end

        @test [(Coefficient(1), [x(1)])] + [(Coefficient(-1), [x(1)])] ≐ 0
    end

    @testset "test multiplication" begin
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test prod(map(algebraElement, nums)) ≐ prod(nums)
        end

        for _ in 1:numRandomTests
            c = randNum()
            a = randAlgebraElement()
            @test c * a ≐ Coefficient(c) * a
        end

        for _ in 1:numRandomTests
            c = randNum()
            a = randAlgebraElement()
            @test c * a ≐ a * c
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test 0 * a ≐ a * 0 ≐ 0
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test 1 * a ≐ a * 1 ≐ a
        end

        for _ in 1:numRandomTests
            c = randNum()
            b = randBasisElement()
            m = randMonomial()
            a = randAlgebraElement()
            @test c * b * a * m ≐ prod(algebraElement, [c, b, a, m])
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            b = randBasisElement()
            @test l * b ≐ [(Coefficient(l), [b])]
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            m = randMonomial()
            @test l * m ≐ [(Coefficient(l), m)]
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            a = randAlgebraElement()
            @test l * a ≐ algebraElement([(Coefficient(l*coeff), mon) for (coeff, mon) in a])
        end
        
        for _ in 1:numRandomTests
            ind = randNums(10)
            @test x(ind) == prod(x, ind)
        end
        
        ind1 = randNums(10)
        ind2 = randNums(10)
        @test x(ind1)*x(ind2) == x([ind1; ind2]) == x(ind1..., ind2...)

        @test x(1) * x(2,3) == x(1,2,3)
        @test !(x(1,2) ≐ x(2,1))
        @test !(x(1) * x(2,3) ≐ x(2,3) * x(1))
        @test (x(1) + x(2)) * (x(1) + x(2)) ≐ x(1,1) + x(1,2) + x(2,1) + x(2,2)
    end

    @testset "test exponentiation" begin
        for _ in 1:numRandomTests
            b = randNum()
            n = randLength(1)
            @test algebraElement(b)^n ≐ b^n
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test a^0 ≐ 1
        end

        for _ in 1:numRandomTests
            n = rand(1:4)
            a = randAlgebraElement()
            @test a^n ≐ prod([a for _ in 1:n])
        end

        for _ in 1:numRandomTests
            n = rand(1:3)
            a = randAlgebraElement()
            @test a*(a^n) ≐ a^(n+1) ≐ (a^n)*a
        end
    end

    @testset "test subtraction" begin
        for _ in 1:numRandomTests
            n = randNum()
            @test -algebraElement(n) == algebraElement(-n)
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test a-a ≐ 0
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test -(-a) ≐ a
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test 0-a ≐ -a
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            b = randAlgebraElement()
            @test -(a+b) ≐ -a-b
        end
        ###
        b = randNum()
        n = randLength(1)
        @test algebraElement(b)^n ≐ b^n

        a = randAlgebraElement()
        @test a^0 ≐ 1

        n = rand(1:4)
        a = randAlgebraElement()
        @test a^n ≐ prod([a for _ in 1:n])

        n = rand(1:3)
        a = randAlgebraElement()
        @test a^(n+1) ≐ a*(a^n)
        @test a^(n+1) ≐ (a^n)*a
    end

    @testset "test subtraction" begin
        n = randNum()
        @test -algebraElement(n) == algebraElement(-n)

        a = randAlgebraElement()
        @test a-a ≐ 0

        a = randAlgebraElement()
        @test -(-a) ≐ a

        a = randAlgebraElement()
        @test 0-a ≐ -a

        a = randAlgebraElement()
        b = randAlgebraElement()
        @test -(a+b) ≐ -a-b
    end

end
