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


function _normalForm(alg::QuadraticAlgebra, ind::Vector{BasisIndex}) :: LinearCombination{Vector{BasisIndex}}
    toIndex(prod::Product{BasisElement}) = map(b -> findfirst(isequal(b), alg.basis), prod) :: Vector{BasisIndex}

    todo = [(1//1, ind)] :: LinearCombination{Vector{BasisIndex}}
    result = LinearCombination{Vector{BasisIndex}}([])

    while !isempty(todo)
        coeff, currInd = pop!(todo)

        changed = false
        for i in 1:length(currInd)-1
            if haskey(alg.relTable, (alg.basis[currInd[i]], alg.basis[currInd[i+1]]))
                linComb = alg.relTable[(alg.basis[currInd[i]], alg.basis[currInd[i+1]])]

                append!(
                    todo,
                    vcat((c * coeff, [currInd[1:i-1]..., toIndex(prod)..., currInd[i+2:end]...]) for (c, prod) in linComb)...,
                )
                changed = true
                break
            end
        end

        if !changed
            push!(result, (coeff, currInd))
        end
    end

    return result
end

function _collect(alg::QuadraticAlgebra, factors::Vector{SymPy.Sym}) :: SymPy.Sym
    factors = [factors...]
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


function normalForm(alg::QuadraticAlgebra, expr::SymPy.Sym) :: SymPy.Sym
    xsum(coll) = isempty(coll) ? toSymPy(0) : sum(coll)

    return expr.
        expand().
        replace(
            sympy.Pow,
            (base, exp) -> _collect(alg, fill(base, fromSymPy(exp)))
        ).
        replace(
            sympy.Mul,
            (m...) -> _collect(alg, [m...])
        ).
        replace(
            f -> f.func == alg.x,
            f -> xsum(
                coeff * alg.x(newInd...)
                for (coeff, newInd) in _normalForm(alg, map(fromSymPy, [f.args...]))
            )
        ).
        simplify()
end


function comm(alg::QuadraticAlgebra, expr1::SymPy.Sym, expr2::SymPy.Sym) :: SymPy.Sym
    return normalForm(alg, expr1 * expr2 - expr2 * expr1)
end

function Base.in(b::BasisElement, alg::QuadraticAlgebra) :: Bool
    return b in alg.basis
end

function Base.in(p::Product{BasisElement}, alg::QuadraticAlgebra) :: Bool
    return all(b -> b in alg, p)
end

function Base.in(a::AlgebraElement, alg::QuadraticAlgebra) :: Bool
    return all(comb -> comb[2] in alg, a)
end

function Base.in(expr::SymPy.Sym, alg::QuadraticAlgebra) :: Bool
    for ex in sympy.preorder_traversal(expr)
        if ex.is_Function && ex.func == alg.x && any(i -> !(i in 1:length(alg.basis)), map(fromSymPy, ex.args))
            return false
        end
    end

    return true
end
