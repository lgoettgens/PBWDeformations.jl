struct SmashProductSymmetricDeformationLie <: AlgebraWithCommutators
    nL :: Int64
    nV :: Int64
    ```
    Contains the simplified commutator of two elements as a linear combination.
    An empty list thus is the empty sum and means that the two elements commutate.
    An absent entry means that there is only a formal commutator.
    ```
    commTable :: Dict{Tuple{BasisElement, BasisElement}, Vector{Tuple{Coefficient, BasisElement}}}

    function SmashProductSymmetricDeformationLie(dynkin::Char, n::Int64, lambda::Vector{Int64})
        @assert n == length(lambda)
        sanitizeLieInput(dynkin, n)

        sp = SmashProductLie(dynkin, n, lambda)
        nL = sp.nL
        nV = sp.nV
        commTable = sp.commTable
        for i in 1:nV, j in 1:i-1
            commTable[(mod(i), mod(j))] = []
        end

        new(nL, nV, commTable)
    end
end
