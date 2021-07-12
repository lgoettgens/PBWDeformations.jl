struct QuadraticAlgebra{T}
    basis :: Vector{BasisElement}

    """
    Stores relations of the form ab = c for basis elements a,b.
    An empty list represents the empty sum, in which case ab = 0.
    An absent entry means that there is no relation, so we cannot simplify ab.
    """
    relTable :: Dict{Tuple{BasisElement, BasisElement}, AlgebraElement}
    extraData :: T
    x :: SymFunction

    QuadraticAlgebra{T}(basis, relTable, extraData = nothing) where T =
        new{T}(basis, relTable, extraData, SymFunction("x", commutative=false))
end

function Base.:(==)(alg1::QuadraticAlgebra, alg2::QuadraticAlgebra) :: Bool
    (alg1.basis, alg1.relTable, alg1.extraData) ==
    (alg2.basis, alg2.relTable, alg2.extraData)
end

function Base.show(io::IO, alg::QuadraticAlgebra) :: Nothing
    println(io, "Algebra with quadratic relations of dimension ", length(alg.basis))
    println(io, "Relation table has ", length(alg.relTable), " entries")
    println(io, "Extra data:")
    print(io, alg.extraData)
end

function Base.in(b::BasisElement, alg::QuadraticAlgebra) :: Bool
    return b in alg.basis
end

function Base.in(m::Monomial{BasisElement}, alg::QuadraticAlgebra) :: Bool
    return all(b -> b in alg, m)
end

function Base.in(a::AlgebraElement, alg::QuadraticAlgebra) :: Bool
    return all(m in alg for (c, m) in a)
end

function Base.in(expr::SymPy.Sym, alg::QuadraticAlgebra) :: Bool
    for ex in sympy.preorder_traversal(expr)
        if ex.is_Function && ex.func == alg.x && any(i -> !(i in 1:length(alg.basis)), map(fromSymPy, ex.args))
            return false
        end
    end

    return true
end


function toIndex(alg::QuadraticAlgebra, b::BasisElement) :: BasisIndex
    return findfirst(isequal(b), alg.basis)
end

function toIndex(alg::QuadraticAlgebra, m::Monomial{BasisElement}) :: Vector{BasisIndex}
    return map(b -> toIndex(alg, b), m)
end


function toSymPy(alg::QuadraticAlgebra, m::Union{BasisElement, Monomial{BasisElement}}) :: SymPy.Sym
    return isempty(m) ? sympify(1) : alg.x(toIndex(alg, m)...)
end

function toSymPy(alg::QuadraticAlgebra, a::AlgebraElement) :: SymPy.Sym
    return sum([c*toSymPy(alg, m) for (c, m) in a], init=sympify(0))
end


function algebraElement(alg::QuadraticAlgebra, i::BasisIndex) :: BasisElement
    return alg.basis[i]
end

function algebraElement(alg::QuadraticAlgebra, ind::Vector{BasisIndex}) :: Monomial{BasisElement}
    return getindex(alg.basis, ind)
end

function algebraElement(alg::QuadraticAlgebra, expr::SymPy.Sym) :: AlgebraElement
    function collectFactors(factors::Vector{SymPy.Sym}) :: SymPy.Sym
        done = false

        while !done
            done = true

            for i in 1:length(factors)-1
                if factors[i].func == alg.x && factors[i+1].func == alg.x
                    done = false

                    factors[i+1] = alg.x(factors[i].args..., factors[i+1].args...)
                    deleteat!(factors, i)
                    break
                end
            end
        end

        return prod(factors)
    end

    function parseAtom(ex::SymPy.Sym) :: AlgebraElement
        coeff = Coefficient(1)

        if ex.func == sympy.Mul
            coeff = Coefficient(fromSymPy(ex.args[1]))
            ex = ex.args[2]
        end

        return coeff * algebraElement(alg, map(fromSymPy, collect(ex.args)))
    end

    @assert expr in alg

    expr = expr.
        expand().
        replace(
            sympy.Pow,
            (base, exp) -> collectFactors(fill(base, fromSymPy(exp)))
        ).
        replace(
            sympy.Mul,
            (m...) -> collectFactors(collect(m))
        )

    if expr.is_number
        return algebraElement(Coefficient(fromSymPy(expr)))
    elseif expr.func in (alg.x, sympy.Mul)
        return parseAtom(expr)
    elseif expr.func == sympy.Add
        return sum(parseAtom, expr.args)
    else
        throw(ArgumentError("cannot parse expr"))
    end
end


function normalForm(alg::QuadraticAlgebra, a::AlgebraElement) :: AlgebraElement
    todo = copy(a)
    result = algebraElement(0)

    while !isempty(todo)
        coeff, mon = pop!(todo)

        changed = false
        for i in 1:length(mon)-1
            if haskey(alg.relTable, (mon[i], mon[i+1]))
                changed = true

                todo += coeff * (mon[1:i-1] * alg.relTable[(mon[i], mon[i+1])] * mon[i+2:end])

                break
            end
        end

        if !changed
            push!(result, (coeff, mon))
        end
    end

    return result
end

function normalForm(alg::QuadraticAlgebra, expr::SymPy.Sym) :: SymPy.Sym
    return toSymPy(alg, normalForm(alg, algebraElement(alg, expr)))
end
