
Coefficient = PD.Coefficient
BasisElement = PD.BasisElement
Monomial{T} = PD.Monomial{T}
AlgebraElement = PD.AlgebraElement

algebraElement = PD.algebraElement

sym(x) = (:sym, x)
coeff(x) = Coefficient(x)
sameSum = PD.sameSum
randNum() = rand(-20:20)
randNums(quantity) = rand(-100:100, quantity)
randLength(start=0) = rand(start:10)
randBasisElement() = sym(randNum()) :: BasisElement
randMonomial() = [randBasisElement() for _ in 1:randLength()] :: Monomial{BasisElement}
randAlgebraElement() = [(randNum()//1, randMonomial()) for _ in 1:randLength()] :: AlgebraElement

@testset ExtendedTestSet "All AlgebraElement.jl tests" begin
    @testset "test conversions" begin
        @test algebraElement() == algebraElement(0) == algebraElement(coeff(0)) == AlgebraElement([]) == AlgebraElement()
        @test algebraElement(algebraElement(0)) == algebraElement(0)

        @test algebraElement(1) == algebraElement(Coefficient(1)) == algebraElement(Monomial{BasisElement}()) == [(coeff(1), Monomial{BasisElement}())] :: AlgebraElement
        @test algebraElement(algebraElement(1)) == algebraElement(1)

        @test algebraElement(sym(1)) == algebraElement([sym(1)]) == [(coeff(1), [sym(1)])] :: AlgebraElement
        @test algebraElement(algebraElement(sym(1))) == algebraElement(sym(1))

        @test algebraElement(map(sym, [1,2,3])) == [(coeff(1), map(sym, [1,2,3]))] :: AlgebraElement
        @test algebraElement(algebraElement(map(sym, [1,2,3]))) == algebraElement(map(sym, [1,2,3]))

        a = [(coeff(1//3), [sym(1)]), (coeff(-1), [sym(1), sym(1)])] :: AlgebraElement
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
        @test sameSum(a, PD.collectSummands(a))

        a = randAlgebraElement()
        @test sameSum(a, shuffle(a))

        c1, c2, c3 = map(coeff, randNums(3))
        mon = randMonomial()
        mon2 = randMonomial()
        @test sameSum([(c1, mon), (c2, mon), (c3, mon2)], [(c3, mon2), (c1+c2, mon)])

        c1, c2 = map(coeff, randNums(2))
        mon = randMonomial()
        mon2 = randMonomial()
        @test sameSum([(c1, mon), (c2, mon2), (-c1, mon)], [(c2, mon2)])
    end

    @testset "test addition" begin
        nums = randNums(randLength(1))
        @test sameSum(sum(map(algebraElement, nums)), algebraElement(sum(nums)))

        as = [randAlgebraElement() for _ in 1:randLength(2)]
        @test sameSum(sum(shuffle(as)), sum(as))

        a = randAlgebraElement()
        @test sameSum(0 + a, a)

        x = randNum()
        a = randAlgebraElement()
        @test sameSum(x + a, coeff(x) + a)

        x = randNum()
        b = randBasisElement()
        m = randMonomial()
        a = randAlgebraElement()
        @test sameSum(x + b + a + m, m + a + x + b)
        @test sameSum(x + b + a + m, sum(map(algebraElement, shuffle([x, b, a, m]))))

        l = randLength(1)
        b = randBasisElement()
        @test sameSum([(coeff(l), [b])], sum(fill(b, l)))

        l = randLength(1)
        m = randMonomial()
        @test sameSum([(coeff(l), m)], sum(fill(m, l)))

        l = randLength(1)
        a = randAlgebraElement()
        @test sameSum([(l*coeff, mon) for (coeff, mon) in a], sum(fill(a, l)))

        @test sameSum(algebraElement(0), [(coeff(1), [sym(1)])] + [(coeff(-1), [sym(1)])])
    end

    @testset "test multiplication" begin
        nums = randNums(randLength(1))
        @test sameSum(prod(map(algebraElement, nums)), algebraElement(prod(nums)))

        c = randNum()
        a = randAlgebraElement()
        @test sameSum(c * a, coeff(c) * a)

        c = randNum()
        a = randAlgebraElement()
        @test sameSum(c * a, a * c)

        a = randAlgebraElement()
        @test sameSum(0 * a, algebraElement(0))

        a = randAlgebraElement()
        @test sameSum(1 * a, a)

        x = randNum()
        b = randBasisElement()
        m = randMonomial()
        a = randAlgebraElement()
        @test sameSum(x * b * a * m, prod(map(algebraElement, [x, b, a, m])))

        l = randLength(1)
        b = randBasisElement()
        @test sameSum(l * b, [(coeff(l), [b])])

        l = randLength(1)
        m = randMonomial()
        @test sameSum(l * m, [(coeff(l), m)])

        l = randLength(1)
        a = randAlgebraElement()
        @test sameSum(l * a, [(l*coeff, mon) for (coeff, mon) in a])

        @test sym(1) * sym(2) != sym(2) * sym(1)
        @test sym(1) * [sym(2), sym(3)] != [sym(2), sym(3)] * sym(1)
        @test sym(1) * [sym(2), sym(3)] != [sym(2), sym(3)] * sym(1)

        @test sameSum((sym(1)+sym(2))*(sym(1)+sym(2)), sym(1)*sym(1) + sym(1)*sym(2) + sym(2)*sym(1) + sym(2)*sym(2))
    end

    @testset "test exponentiation" begin
        b = randNum()
        n = randLength(1)
        @test sameSum(algebraElement(b)^n, algebraElement(b^n))

        a = randAlgebraElement()
        @test sameSum(a^0, algebraElement(1))

        n = rand(1:3)
        a = randAlgebraElement()
        @test sameSum(a^n, prod([a for _ in 1:n]))

        n = rand(1:3)
        a = randAlgebraElement()
        @test sameSum(a^(n+1), a*(a^n))
    end

    @testset "test subtraction" begin
        n = randNum()
        @test sameSum(-algebraElement(n), algebraElement(-n))

        a = randAlgebraElement()
        @test sameSum(a-a, algebraElement(0))

        a = randAlgebraElement()
        @test sameSum(-(-a), a)

        a = randAlgebraElement()
        @test sameSum(algebraElement(0)-a, -a)

        a = randAlgebraElement()
        b = randAlgebraElement()
        @test sameSum(-(a+b), -a-b)
    end

end
