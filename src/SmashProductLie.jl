struct SmashProductLie
    dynkin :: Char
    n :: Int64
    lambda :: Vector{Int64}
    nL :: Int64
    nV :: Int64
end

function smashProductLie(dynkin::Char, n::Int64, lambda::Vector{Int64}) :: AlgebraWithCommutators{SmashProductLie}
    @assert n == length(lambda)
    sanitizeLieInput(dynkin, n)

    commTable = Dict{Tuple{BasisElement, BasisElement}, LinearCombination}()
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

    extraData = SmashProductLie(dynkin, n, lambda, nL, nV)
    basis = [[lie(i) for i in 1:nL]..., [mod(i) for i in 1:nV]...] :: Vector{BasisElement}
    return AlgebraWithCommutators{SmashProductLie}(basis, commTable, extraData)
end
