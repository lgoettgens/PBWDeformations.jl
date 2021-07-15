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
        nums = randNums(randLength(1))
        @test PD.collectSummands(map(algebraElement, nums)) == algebraElement(sum(nums))

        b = randBasisElement()
        @test PD.collectSummands(algebraElement(b)) == algebraElement(b)

        mon = randMonomial()
        @test PD.collectSummands(algebraElement(mon)) == algebraElement(mon)

        a = randAlgebraElement()
        @test a ≐ PD.collectSummands(a)

        a = randAlgebraElement()
        @test a ≐ shuffle(a)

        c1, c2, c3 = map(Coefficient, randNums(3))
        mon = randMonomial()
        mon2 = randMonomial()
        @test [(c1, mon), (c2, mon), (c3, mon2)] ≐ [(c3, mon2), (c1+c2, mon)]

        c1, c2 = map(Coefficient, randNums(2))
        mon = randMonomial()
        mon2 = randMonomial()
        @test [(c1, mon), (c2, mon2), (-c1, mon)] ≐ [(c2, mon2)]
    end

    @testset "test addition" begin
        nums = randNums(randLength(1))
        @test sum(map(algebraElement, nums)) ≐ algebraElement(sum(nums))

        as = [randAlgebraElement() for _ in 1:randLength(2)]
        @test sum(shuffle(as)) ≐ sum(as)

        a = randAlgebraElement()
        @test 0 + a ≐ a

        c = randNum()
        a = randAlgebraElement()
        @test c + a ≐ Coefficient(c) + a

        c = randNum()
        b = randBasisElement()
        m = randMonomial()
        a = randAlgebraElement()
        @test c + b + a + m ≐ sum(algebraElement, shuffle([c, b, a, m]))

        l = randLength(1)
        b = randBasisElement()
        @test [(Coefficient(l), [b])] ≐ sum(fill(b, l))

        l = randLength(1)
        m = randMonomial()
        @test [(Coefficient(l), m)] ≐ sum(fill(m, l))

        l = randLength(1)
        a = randAlgebraElement()
        b = [(Coefficient(l*coeff), mon) for (coeff, mon) in a] :: AlgebraElement
        @test sum(fill(a,l)) ≐ b

        @test [(Coefficient(1), [x(1)])] + [(Coefficient(-1), [x(1)])] ≐ 0
    end

    @testset "test multiplication" begin
        nums = randNums(randLength(1))
        @test prod(map(algebraElement, nums)) ≐ prod(nums)

        c = randNum()
        a = randAlgebraElement()
        @test c * a ≐ Coefficient(c) * a

        c = randNum()
        a = randAlgebraElement()
        @test c * a ≐ a * c

        a = randAlgebraElement()
        @test 0 * a ≐ 0

        a = randAlgebraElement()
        @test 1 * a ≐ a

        c = randNum()
        b = randBasisElement()
        m = randMonomial()
        a = randAlgebraElement()
        @test c * b * a * m ≐ prod(algebraElement, [c, b, a, m])

        l = randLength(1)
        b = randBasisElement()
        @test l * b ≐ [(Coefficient(l), [b])]

        l = randLength(1)
        m = randMonomial()
        @test l * m ≐ [(Coefficient(l), m)]

        l = randLength(1)
        a = randAlgebraElement()
        @test l * a ≐ [(Coefficient(l*coeff), mon) for (coeff, mon) in a]

        ind = randNums(10)
        @test x(ind) == prod(x, ind)

        ind1 = randNums(10)
        ind2 = randNums(10)
        @test x(ind1)*x(ind2) == x([ind1; ind2]) == x(ind1..., ind2...)

        @test x(1) * x(2,3) == x(1,2,3)
        @test !(x(1,2) ≐ x(2,1))
        @test !(x(1) * x(2,3) ≐ x(2,3) * x(1))
        @test (x(1) + x(2)) * (x(1) + x(2)) ≐ x(1,1) + x(1,2) + x(2,1) + x(2,2)
    end

    @testset "test exponentiation" begin
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
