struct QuadraticAlgebra{T}
    basis :: Vector{BasisElement}

    """
    #FIXME
    An empty list thus is the empty sum and means that the two elements commute.
    An absent entry means that there is only a formal commutator.
    """
    relTable :: Dict{Tuple{BasisElement, BasisElement}, AlgebraElement}
    extraData :: T
    x :: SymFunction

    QuadraticAlgebra{T}(basis, relTable, extraData = nothing) where T = new{T}(basis, relTable, extraData, SymFunction("x"))
end

function Base.:(==)(alg1::QuadraticAlgebra, alg2::QuadraticAlgebra)
    # also check for x?
    (alg1.basis, alg1.relTable, alg1.extraData) ==
    (alg2.basis, alg2.relTable, alg2.extraData)
end

function Base.show(io::IO, alg::QuadraticAlgebra)
    println(io, "Algebra with quadratic relations of dimension ", length(alg.basis))
    println(io, "Relation table has ", length(alg.relTable), " entries")
    println(io, "Extra data:")
    show(alg.extraData)
end


function _normalForm(alg::QuadraticAlgebra, ind::Vector{BasisIndex}) :: LinearCombination{Vector{BasisIndex}}
    toIndex(prod :: Product{BasisElement}) = map(b -> findfirst(isequal(b), alg.basis), prod) :: Vector{BasisIndex}

    todo = [(1, ind)] :: LinearCombination{Vector{BasisIndex}}
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

function normalForm(alg::QuadraticAlgebra, expr::SymPy.Sym) :: SymPy.Sym
    xsum(coll) = isempty(coll) ? 0 : sum(coll)

    expr.replace(
        f -> f.func == alg.x,
        f -> xsum(
            coeff * alg.x(newInd...)
            for (coeff, newInd) in _normalForm(alg, map(fromSymPy, [f.args...]))
        )
    )
end
