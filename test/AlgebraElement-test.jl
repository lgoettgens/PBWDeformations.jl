sym(x) = (:sym, x)
coeff(x) = Coefficient(x)
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
            @test sameSum(a, PD.collectSummands(a))
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test sameSum(a, shuffle(a))
        end

        for _ in 1:numRandomTests
            c1, c2, c3 = map(coeff, randNums(3))
            mon = randMonomial()
            mon2 = randMonomial()
            @test sameSum([(c1, mon), (c2, mon), (c3, mon2)], [(c3, mon2), (c1+c2, mon)])
        end

        for _ in 1:numRandomTests
            c1, c2 = map(coeff, randNums(2))
            mon = randMonomial()
            mon2 = randMonomial()
            @test sameSum([(c1, mon), (c2, mon2), (-c1, mon)], [(c2, mon2)])
        end
    end

    @testset "test addition" begin
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test sameSum(sum(map(algebraElement, nums)), algebraElement(sum(nums)))
        end

        for _ in 1:numRandomTests
            as = [randAlgebraElement() for _ in 1:randLength(2)]
            @test sameSum(sum(shuffle(as)), sum(as))
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test sameSum(0 + a, a)
        end

        for _ in 1:numRandomTests
            x = randNum()
            a = randAlgebraElement()
            @test sameSum(x + a, coeff(x) + a)
        end

        for _ in 1:numRandomTests
            x = randNum()
            b = randBasisElement()
            m = randMonomial()
            a = randAlgebraElement()
            @test sameSum(x + b + a + m, m + a + x + b)
            @test sameSum(x + b + a + m, sum(map(algebraElement, shuffle([x, b, a, m]))))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            b = randBasisElement()
            @test sameSum([(coeff(l), [b])], sum(fill(b, l)))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            m = randMonomial()
            @test sameSum([(coeff(l), m)], sum(fill(m, l)))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            a = randAlgebraElement()
            b = algebraElement([(Coefficient(l*coeff), mon) for (coeff, mon) in a])
            @test sameSum(sum(fill(a,l)), b)
        end

        @test sameSum(algebraElement(0), [(coeff(1), [sym(1)])] + [(coeff(-1), [sym(1)])])
    end

    @testset "test multiplication" begin
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test sameSum(prod(map(algebraElement, nums)), algebraElement(prod(nums)))
        end

        for _ in 1:numRandomTests
            c = randNum()
            a = randAlgebraElement()
            @test sameSum(c * a, coeff(c) * a)
        end

        for _ in 1:numRandomTests
            c = randNum()
            a = randAlgebraElement()
            @test sameSum(c * a, a * c)
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test sameSum(0 * a, algebraElement(0))
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test sameSum(1 * a, a)
        end

        for _ in 1:numRandomTests
            x = randNum()
            b = randBasisElement()
            m = randMonomial()
            a = randAlgebraElement()
            @test sameSum(x * b * a * m, prod(map(algebraElement, [x, b, a, m])))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            b = randBasisElement()
            @test sameSum(l * b, [(coeff(l), [b])])
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            m = randMonomial()
            @test sameSum(l * m, [(coeff(l), m)])
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            a = randAlgebraElement()
            @test sameSum(l * a, algebraElement([(Coefficient(l*coeff), mon) for (coeff, mon) in a]))
        end

        @test sym(1) * sym(2) != sym(2) * sym(1)
        @test sym(1) * [sym(2), sym(3)] != [sym(2), sym(3)] * sym(1)
        @test sym(1) * [sym(2), sym(3)] != [sym(2), sym(3)] * sym(1)

        @test sameSum((sym(1)+sym(2))*(sym(1)+sym(2)), sym(1)*sym(1) + sym(1)*sym(2) + sym(2)*sym(1) + sym(2)*sym(2))
    end

    @testset "test exponentiation" begin
        for _ in 1:numRandomTests
            b = randNum()
            n = randLength(1)
            @test sameSum(algebraElement(b)^n, algebraElement(b^n))
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test sameSum(a^0, algebraElement(1))
        end

        for _ in 1:numRandomTests
            n = rand(1:3)
            a = randAlgebraElement()
            @test sameSum(a^n, prod([a for _ in 1:n]))
        end

        for _ in 1:numRandomTests
            n = rand(1:3)
            a = randAlgebraElement()
            @test sameSum(a^(n+1), a*(a^n))
        end
    end

    @testset "test subtraction" begin
        for _ in 1:numRandomTests
            n = randNum()
            @test sameSum(-algebraElement(n), algebraElement(-n))
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test sameSum(a-a, algebraElement(0))
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test sameSum(-(-a), a)
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            @test sameSum(algebraElement(0)-a, -a)
        end

        for _ in 1:numRandomTests
            a = randAlgebraElement()
            b = randAlgebraElement()
            @test sameSum(-(a+b), -a-b)
        end
    end

end
