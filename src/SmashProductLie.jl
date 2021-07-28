struct SmashProductLie
    dynkin :: Char
    n :: Int64
    lambda :: Vector{Int64}
    nL :: Int64
    nV :: Int64
    matrixRepL :: Vector{Matrix{Int64}}
end

function Base.:(==)(sp1::SmashProductLie, sp2::SmashProductLie) :: Bool
    (sp1.dynkin, sp1.n, sp1.lambda, sp1.nL, sp1.nV, sp1.matrixRepL) ==
    (sp2.dynkin, sp2.n, sp2.lambda, sp2.nL, sp2.nV, sp2.matrixRepL)
end

function Base.show(io::IO, sp::SmashProductLie) :: Nothing
    println(io, "Smash product of lie algebra with highest weight module")
    println(io, "Lie algebra: type ", sp.dynkin, sp.n, ", dimension ", sp.nL)
    println(io, "Module: highest weight ", sp.lambda, ", dimension ", sp.nV)
end


function smashProductLie(dynkin::Char, n::Int64, lambda::Vector{Int64}, C::Type = Rational{Int64}) :: QuadraticAlgebra{C, SmashProductLie}
    @assert n == length(lambda)
    sanitizeLieInput(dynkin, n)

    relTable = Dict{Tuple{BasisElementInternal, BasisElementInternal}, AlgebraElementInternal{C}}()
        # (lie(i), lie(j)) => [(c, [lie(k)])]
        # (lie(i), mod(j)) => [(c, [mod(k)])]

    L = GAP.SimpleLieAlgebra(toGAP(string(dynkin)), n, GAP.Rationals)
    nL = GAP.Dimension(L)
    bL = GAP.BasisVectors(GAP.Basis(L))
    commTableL = fromGAP(GAP.StructureConstantsTable(GAP.Basis(L)))[1:nL]
    matrixRepL = getMatrixRep(L, [i == 1 ? 1 : 0 for i in 1:n])

    for i in 1:nL, j in 1:i-1
        relTable[(lieInt(i), lieInt(j))] = [
            (1, [lieInt(j), lieInt(i)]);
            collect((C(c), [lieInt(k)]) for (k, c) in zip(commTableL[i][j]...));
        ]
    end

    V = GAP.HighestWeightModule(L, toGAP(lambda))
    nV = GAP.Dimension(V)
    bV = GAP.BasisVectors(GAP.Basis(V))

    for i in 1:nL, j in 1:nV
        relTable[(lieInt(i), modInt(j))] = [
            (1, [modInt(j), lieInt(i)]);
            collect((C(c), [modInt(k)]) for (k, c) in enumerate(fromGAP(GAP.Coefficients(GAP.Basis(V), bL[i]^bV[j]))) if !iszero(c));
        ]
    end

    extraData = SmashProductLie(dynkin, n, lambda, nL, nV, matrixRepL)
    basis = [[mod(i) for i in 1:nV]; [lie(i) for i in 1:nL]] :: Vector{BasisElement}

    return QuadraticAlgebra{C, SmashProductLie}(basis, relTable, extraData)
end


function getMatrixRep(dynkin::Char, n::Int64) :: Vector{Matrix{Int64}}
    sanitizeLieInput(dynkin, n)
    
    L = GAP.SimpleLieAlgebra(toGAP(string(dynkin)), n, GAP.Rationals)
    lambda = [i == 1 ? 1 : 0 for i in 1:n]

    return getMatrixRep(L, lambda)
end

function getMatrixRep(L #= :: Gap.LieAlgebra =#, lambda::Vector{Int64}) :: Vector{Matrix{Int64}}
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
