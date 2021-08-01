x = test
randNum() = rand(-20:20)
randNums(quantity) = rand(-20:20, quantity)
randLength(start=0) = rand(start:10)
randBasisElement() = x(randNum()) :: BasisElement{Rational{Int64}}
basis = [x(i) for i in -20:20]
randMonomial() = randMonomial(basis)
randAlgebraElement() = randAlgebraElement(basis)
#randAlgebraElement() = AlgebraElement{Rational{Int64}}()

@testset ExtendedTestSet "All AlgebraElement.jl tests" begin
    @testset "test conversions, Coeff=$C" for C in [Rational{Int64}]
        @test AlgebraElement{C}() == AlgebraElement{C}(0) == AlgebraElement{C}(C(0))
        @test AlgebraElement{C}(AlgebraElement{C}(0)) == AlgebraElement{C}(0)

        @test AlgebraElement{C}(1) == AlgebraElement{C}(Monomial{C}()) == AlgebraElement{C}([(C(1), Monomial{C}())])
        @test AlgebraElement{C}(AlgebraElement{C}(1)) == AlgebraElement{C}(1)

        @test AlgebraElement{C}(x(1; C)) == AlgebraElement{C}(Monomial{C}(x(1; C))) == AlgebraElement{C}(x([1]; C)) == AlgebraElement{C}([(C(1), x([1]; C))])
        @test AlgebraElement{C}(AlgebraElement{C}(x(1; C))) == AlgebraElement{C}(x(1; C))

        @test AlgebraElement{C}(x(1,2,3; C)) == AlgebraElement{C}([(C(1), x(1,2,3; C))])
        @test AlgebraElement{C}(AlgebraElement{C}(x(1,2,3; C))) == AlgebraElement{C}(x(1,2,3; C))

        a = AlgebraElement{C}([(C(1//3), x([1]; C)), (C(-1), x(1,1; C))])
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
            @test a ≐ AlgebraElement{C}(shuffle(unpack(a))) # TODO make shuffle work on wrapper
        end

        for _ in 1:numRandomTests
            c1, c2, c3 = map(C, randNums(3))
            mon = randMonomial()
            mon2 = randMonomial()
            @test AlgebraElement{C}([(c1, mon), (c2, mon), (c3, mon2)]) ≐ AlgebraElement{C}([(c3, mon2), (c1+c2, mon)])
        end

        for _ in 1:numRandomTests
            c1, c2 = map(C, randNums(2))
            mon = randMonomial()
            mon2 = randMonomial()
            @test AlgebraElement{C}([(c1, mon), (c2, mon2), (-c1, mon)]) ≐ AlgebraElement{C}([(c2, mon2)])
        end
    end

    @testset "test addition, Coeff=$C" for C in [Rational{Int64}]
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test sum(map(AlgebraElement{C}, nums)) ≐ AlgebraElement{C}(sum(nums))
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
            @test c + b + a + m ≐ sum(AlgebraElement{C}, shuffle([c, b, a, m]))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            b = randBasisElement()
            @test AlgebraElement{C}([(C(l), Monomial{C}(b))]) ≐ sum(fill(b, l))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            m = randMonomial()
            @test AlgebraElement{C}([(C(l), m)]) ≐ sum(fill(m, l))
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            a1 = randAlgebraElement()
            a2 = AlgebraElement{C}([(C(l*coeff), mon) for (coeff, mon) in a1])
            @test sum(fill(a1,l)) ≐ a2
        end

        @test AlgebraElement{C}([(C(1), x([1]; C))]) + AlgebraElement{C}([(C(-1), x([1]; C))]) ≐ 0
    end

    @testset "test multiplication, Coeff=$C" for C in [Rational{Int64}]
        for _ in 1:numRandomTests
            nums = randNums(randLength(1))
            @test prod(map(AlgebraElement{C}, nums)) ≐ prod(nums)
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
            @test c * b * a * m ≐ prod(AlgebraElement{C}, [c, b, a, m])
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            b = randBasisElement()
            @test l * b ≐ AlgebraElement{C}([(C(l), Monomial{C}(b))])
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            m = randMonomial()
            @test l * m ≐ AlgebraElement{C}([(C(l), m)])
        end

        for _ in 1:numRandomTests
            l = randLength(1)
            a = randAlgebraElement()
            @test l * a ≐ AlgebraElement{C}([(C(l*coeff), mon) for (coeff, mon) in a])
        end
        
        for _ in 1:numRandomTests
            ind = randNums(10)
            @test x(ind; C) == prod(i -> x(i; C), ind)
        end
        
        ind1 = randNums(10)
        ind2 = randNums(10)
        @test x(ind1; C) * x(ind2; C) == x([ind1; ind2]; C) == x(ind1..., ind2...; C)

        @test x([1]; C) == Monomial{C}(x(1; C))
        @test x(1; C) * x(2; C) == x([1, 2]; C) == x(1,2; C) == Monomial{C}(x(1; C), x(2; C)) == Monomial{C}([x(1; C), x(2; C)])
        @test x(1; C) * x(2,3; C) == x(1,2,3; C)
        @test !(x(1,2; C) ≐ x(2,1; C))
        @test !(x(1; C) * x(2,3; C) ≐ x(2,3; C) * x(1; C))
        @test (x(1; C) + x(2; C)) * (x(1; C) + x(2; C)) ≐ x(1,1; C) + x(1,2; C) + x(2,1; C) + x(2,2; C)
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
            n = rand(0:2) # TODO: increase performance so that we have 1:3 here again, as before
            a = randAlgebraElement()
            @test a^n ≐ prod([a for _ in 1:n]; init=AlgebraElement{C}(1))
        end

        for _ in 1:numRandomTests
            n = rand(0:2) # TODO: increase performance so that we have 1:3 here again, as before
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
