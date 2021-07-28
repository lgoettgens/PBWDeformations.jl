x = test
xInt = PD.testInt
randNum() = rand(-20:20)
randNums(quantity) = rand(-20:20, quantity)
randLength(start=0) = rand(start:10)
randBasisElement() = x(randNum()) :: BasisElement
basis = x(-20:20)
randMonomial() = randMonomial(basis)
randAlgebraElement() = randAlgebraElement(basis)

@testset ExtendedTestSet "All AlgebraElement.jl tests" begin
    @testset "test conversions, Coeff=$C" for C in [Rational{Int64}]
        @test AlgebraElement{C}() == AlgebraElement{C}(0) == AlgebraElement{C}(C(0)) == AlgebraElement{C}(PD.AlgebraElementInternal{C}([]))
        @test AlgebraElement{C}(AlgebraElement{C}(0)) == AlgebraElement{C}(0)

        @test AlgebraElement{C}(1) == AlgebraElement{C}(Monomial()) == AlgebraElement{C}([(C(1), Monomial())])
        @test AlgebraElement{C}(AlgebraElement{C}(1)) == AlgebraElement{C}(1)

        @test AlgebraElement{C}(xInt(1)) == AlgebraElement{C}(x(1)) == AlgebraElement{C}([x(1)]) == AlgebraElement{C}([(C(1), [x(1)])])
        @test AlgebraElement{C}(AlgebraElement{C}(x(1))) == AlgebraElement{C}(x(1))

        @test AlgebraElement{C}(x(1,2,3)) == AlgebraElement{C}([(C(1), x(1,2,3))])
        @test AlgebraElement{C}(AlgebraElement{C}(x(1,2,3))) == AlgebraElement{C}(x(1,2,3))

        a = AlgebraElement{C}([(C(1//3), [x(1)]), (C(-1), x(1,1))])
        @test AlgebraElement{C}(a) == a
    end

    @testset "test misc operations, Coeff=$C" for C in [Rational{Int64}]
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test PD.collectSummands(map(AlgebraElement{C}, nums)) == AlgebraElement{C}(sum(nums))
        end

        for _ in 1:numRandomTests
            b = randBasisElement()
            @test PD.collectSummands(AlgebraElement{C}(b)) == AlgebraElement{C}(b)
        end

        for _ in 1:numRandomTests
            mon = randMonomial()
            @test PD.collectSummands(AlgebraElement{C}(mon)) == AlgebraElement{C}(mon)
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
            c1, c2, c3 = map(C, randNums(3))
            mon = randMonomial()
            mon2 = randMonomial()
            @test [(c1, mon), (c2, mon), (c3, mon2)] ≐ [(c3, mon2), (c1+c2, mon)]
        end

        for _ in 1:numRandomTests
            c1, c2 = map(C, randNums(2))
            mon = randMonomial()
            mon2 = randMonomial()
            @test [(c1, mon), (c2, mon2), (-c1, mon)] ≐ [(c2, mon2)]
        end
    end

    @testset "test addition, Coeff=$C" for C in [Rational{Int64}]
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test sum(map(z ->AlgebraElement{C}(z), nums)) ≐ AlgebraElement{C}(sum(nums))
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
            @test c + a ≐ C(c) + a
        end

        for _ in 1:numRandomTests
            c = randNum()
            b = randBasisElement()
            m = randMonomial()
            a = randAlgebraElement()
            @test c + b + a + m ≐ sum(z -> AlgebraElement{C}(z), shuffle([c, b, a, m]))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            b = randBasisElement()
            @test [(C(l), [b])] ≐ sum(fill(b, l))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            m = randMonomial()
            @test [(C(l), m)] ≐ sum(fill(m, l))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            a = randAlgebraElement()
            b = AlgebraElement{C}([(C(l*coeff), mon) for (coeff, mon) in a])
            @test sameSum(sum(fill(a,l)), b)
        end

        @test [(C(1), [x(1)])] + [(C(-1), [x(1)])] ≐ 0
    end

    @testset "test multiplication, Coeff=$C" for C in [Rational{Int64}]
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test prod(map(z -> AlgebraElement{C}(z), nums)) ≐ prod(nums)
        end

        for _ in 1:numRandomTests
            c = randNum()
            a = randAlgebraElement()
            @test c * a ≐ C(c) * a
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
            @test c * b * a * m ≐ prod(z -> AlgebraElement{C}(z), [c, b, a, m])
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            b = randBasisElement()
            @test l * b ≐ [(C(l), [b])]
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            m = randMonomial()
            @test l * m ≐ [(C(l), m)]
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            a = randAlgebraElement()
            @test l * a ≐ AlgebraElement{C}([(C(l*coeff), mon) for (coeff, mon) in a])
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

    @testset "test exponentiation, Coeff=$C" for C in [Rational{Int64}]
        for _ in 1:numRandomTests
            b = randNum()
            n = randLength(1)
            @test AlgebraElement{C}(b)^n ≐ b^n
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

    @testset "test subtraction, Coeff=$C" for C in [Rational{Int64}]
        for _ in 1:numRandomTests
            n = randNum()
            @test -AlgebraElement{C}(n) == AlgebraElement{C}(-n)
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
            @test a-0 ≐ a
        end

        for _ in 1:numRandomTests
            c = randNum()
            a = randAlgebraElement()
            @test a-c ≐ -(c-a)
        end

        for _ in 1:numRandomTests
            b = randBasisElement()
            a = randAlgebraElement()
            @test a-b ≐ -(b-a)
        end

        for _ in 1:numRandomTests
            m = randMonomial()
            a = randAlgebraElement()
            @test a-m ≐ -(m-a)
        end

        for _ in 1:numRandomTests
            a1 = randAlgebraElement()
            a2 = randAlgebraElement()
            @test a1-a2 ≐ -(a2-a1)
        end

        for _ in 1:numRandomTests
            a1 = randAlgebraElement()
            a2 = randAlgebraElement()
            @test -(a1+a2) ≐ -a1-a2
        end
    end

end
