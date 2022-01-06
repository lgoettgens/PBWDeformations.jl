struct SmashProductLie
    dynkin :: Char
    n :: Int64
    lambda :: Vector{Int64}
    nL :: Int64
    nV :: Int64
    matrixRepL :: Vector{Matrix{Int64}}
end

function get_matrix_rep(dynkin::Char, n::Int64) :: Vector{Matrix{Int64}}
    sanitize_lie_input(dynkin, n)
    
    L = GAP.SimpleLieAlgebra(toGAP(string(dynkin)), n, GAP.Rationals)
    lambda = [i == 1 ? 1 : 0 for i in 1:n]

    return get_matrix_rep(L, lambda)
end

function get_matrix_rep(L #= :: Gap.LieAlgebra =#, lambda::Vector{Int64}) :: Vector{Matrix{Int64}}
    @assert GAP.IsLieAlgebra(L)

    nL = GAP.Dimension(L)
    bL = GAP.BasisVectors(GAP.Basis(L))

    V = GAP.HighestWeightModule(L, toGAP(lambda))
    nV = GAP.Dimension(V)
    bV = GAP.BasisVectors(GAP.Basis(V))

    matrixRep = [zeros(Int64, nV,nV) for _ in 1:nL]
    for i in 1:nL, j in 1:nV
        matrixRep[i][:,j] = fromGAP(GAP.Coefficients(GAP.Basis(V), bL[i]^bV[j]))
    end

    return matrixRep
end
