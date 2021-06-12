struct AlgebraWithCommutators{T}
    basis :: Vector{BasisElement}

    """
    Contains the simplified commutator of two elements as a linear combination.
    An empty list thus is the empty sum and means that the two elements commutate.
    An absent entry means that there is only a formal commutator.
    """
    commTable :: Dict{Tuple{BasisElement, BasisElement}, LinearCombination}
    extraData :: T
end

x = sympy.Function("x")

function _normalForm(alg::AlgebraWithCommutators, coeff::Coefficient, ind::Vector{Int64}) :: Vector{Tuple{Coefficient,Vector{Int64}}}
    for i in 1:length(ind)-1
        if haskey(alg.commTable, (alg.basis[ind[i]], alg.basis[ind[i+1]]))
            comm = alg.commTable[(alg.basis[ind[i]], alg.basis[ind[i+1]])]

            return [
                _normalForm(alg, coeff, [ind[1:i-1]..., ind[i+1], ind[i], ind[i+2:end]...])...,
                vcat(
                    (_normalForm(alg, c * coeff, [ind[1:i-1]..., findfirst(isequal(b), alg.basis), ind[i+2:end]...]) for (c, b) in comm)...
                )...
            ]            
        end
    end
    return [(coeff, ind)]
end


function normalForm(alg, expr)
    xsum(coll) = isempty(coll) ? 0 : sum(coll)

    expr.replace(
        f -> f.func == x, 
        f -> xsum(
            coeff * x(newInd...)
            for (coeff, newInd) in _normalForm(alg, 1, map(fromSymPy, [f.args...]))
        )
    )
end
