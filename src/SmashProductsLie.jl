using Oscar

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

function sanitizeLieInput(dynkin::Char, n::Integer) :: Nothing
    @assert dynkin in ['A', 'B', 'C', 'D']
    if dynkin == 'A'
        @assert n >= 1
    elseif dynkin == 'B' || dykin == 'C'
        @assert n >= 2
    elseif dynkin == 'D'
        @assert n >= 4
    end
end

BasisElement = Pair{Symbol, Integer}
Coefficient = Integer

lie(i::Integer) = (:lie, i) :: BasisElement
mod(i::Integer) = (:mod, i) :: BasisElement

struct SmashProductLie <: SmashProduct
    nL :: Integer
    nV :: Integer
    ```
    Contains the simplified commutator of two elements as a linear combination.
    An empty list thus is the empty sum and means that the two elements commutate.
    An absent entry means that there is only a formal commutator.
    ```
    commTable :: Dict{Pair{BasisElement, BasisElement}, Vector{Pair{Coefficient, BasisElement}}}

    function SmashProductLie(dynkin::Char, n::Integer, lambda::Vector{Integer})
        @assert n == length(lambda)
        sanitizeLieInput(dynkin, n)

        commTable = Dict{Pair{BasisElement, BasisElement}, Vector{Pair{Coefficient, BasisElement}}}()
            # (lie(i), lie(j)) => [(c, lie(k))]
            # (lie(i), mod(j)) => [(c, mod(k))]

        L = GAP.SimpleLieAlgebra(toGAP(string(dynkin)), n, GAP.Rationals)
        nL = GAP.Dimension(L)
        commTableL = fromGAP(GAP.StructureConstantsTable(GAP.Basis(L)))[1:nL]
        
        for i in 1:nL, j in 1:i-1
            commTable[(lie(i), lie(j))] = [(c, lie(k)) for (k, c) in zip(commTableL[i][j]...)]
        end

        V = GAP.HighestWeightModule(L, toGAP(lambda))
        nV = GAP.Dimension(V)
        bL = GAP.BasisVectors(GAP.Basis(L))
        bV = GAP.BasisVectors(GAP.Basis(V))

        for i in 1:nL, j in 1:nV
            commTable[(lie(i), mod(j))] = [
                (c, mod(k)) 
                for (k, c) in enumerate(fromGAP(GAP.Coefficients(GAP.Basis(V), bL[i]^bV[j])))
                if c != 0
            ]
        end

        new(nL, nV, commTable)
    end
end
