struct SmashProductDeformLie{C}
    sp :: SmashProductLie
    symmetric :: Bool
    kappa :: Matrix{AlgebraElement{C}}
end

function Base.:(==)(spd1::SmashProductDeformLie{C}, spd2::SmashProductDeformLie{C}) :: Bool where C
    (spd1.sp, spd1.symmetric) == (spd2.sp, spd2.symmetric)
end

function Base.show(io::IO, spd::SmashProductDeformLie) :: Nothing
    if spd.symmetric
        println(io, "Symmetric deformation of:")
    else
        println(io, "Deformation of:")
    end
    print(io, spd.sp)
end


function smashProductDeformLie(sp::QuadraticAlgebra{C, SmashProductLie}, kappa::Matrix{AlgebraElement{C}}, one=C(1)::C) :: QuadraticAlgebra{C, SmashProductDeformLie{C}} where C
    nV = sp.extraData.nV
    @assert size(kappa) == (nV, nV) "size of kappa matches module dimension"

    # basis of smash product consists of basis of module and basis of Hopf algebra
    hopfBasis = filter(!ismod, sp.basis)
    @assert all(e -> issubset(basisElements(e), hopfBasis), kappa) "kappa only takes values in Hopf algebra"

    for i in 1:nV, j in 1:i
        @assert kappa[i,j] ≐ -kappa[j,i] "kappa is skew-symmetric"
    end

    relTable = sp.relTable
    symmetric = true

    for i in 1:nV, j in 1:i-1
        symmetric &= (kappa[i,j] ≐ 0)

        # We have the commutator relation [mod(i), mod(j)] = kappa[i,j]
        # which is equivalent to mod(i)*mod(j) = mod(j)*mod(i) + kappa[i,j]
        relTable[(mod(i; C), mod(j; C))] = AlgebraElement{C}(mod(j, i; C), one) + kappa[i,j]
    end

    extraData = SmashProductDeformLie{C}(sp.extraData, symmetric, kappa)
    return QuadraticAlgebra{C, SmashProductDeformLie{C}}(sp.basis, relTable, extraData)
end


function smashProductSymmDeformLie(sp::QuadraticAlgebra{C, SmashProductLie}) :: QuadraticAlgebra{C, SmashProductDeformLie{C}} where C
    relTable = sp.relTable

    for i in 1:sp.extraData.nV, j in 1:i-1
        relTable[(mod(i; C), mod(j; C))] = AlgebraElement{C}(mod(j, i; C))
    end

    extraData = SmashProductDeformLie{C}(sp.extraData, true, fill(AlgebraElement{C}(), sp.extraData.nV, sp.extraData.nV))
    return QuadraticAlgebra{C, SmashProductDeformLie{C}}(sp.basis, relTable, extraData)
end

function smashProductSymmDeformLie(dynkin::Char, n::Int64, lambda::Vector{Int64}; C::Type=Rational{Int64}) :: QuadraticAlgebra{C, SmashProductDeformLie{C}}
    @assert n == length(lambda)
    sanitizeLieInput(dynkin, n)

    return smashProductSymmDeformLie(smashProductLie(dynkin, n, lambda; C))
end

struct PBWDeformEqs{C}
    d :: QuadraticAlgebra{C, SmashProductDeformLie{C}}
    one :: C

    PBWDeformEqs{C}(d::QuadraticAlgebra{C, SmashProductDeformLie{C}}, one = C(1)) where C = new{C}(d, one)
end

function Base.iterate(eqs::PBWDeformEqs{C}, s=nothing) where C
    nL = eqs.d.extraData.sp.nL
    nV = eqs.d.extraData.sp.nV
    kappa = eqs.d.extraData.kappa

    #s = (flag = 1, c_state, h) oder (flag = 2, i, j, k)

    if s === nothing # may fail for nV<2
        s = Any[1, 1, nothing]
    end

    if s[1] == 1
        comb = Combinatorics.Combinations(nV,2)
        res = s[3] === nothing ? iterate(comb) : iterate(comb, s[3])
        if res === nothing
            s[2] += 1
            if s[2] > nL
                s[1] = 2
                s[3] = nothing
            else
                comb = Combinatorics.Combinations(nV,2)
                res = iterate(comb)
                s[3] = res[2]
            end
        else
            s[3] = res[2]
        end
    end
    if s[1] == 2
        comb = Combinatorics.Combinations(nV,3)
        res = s[3] === nothing ? iterate(comb) : iterate(comb, s[3])
        if res === nothing
            return nothing
        else
            s[3] = res[2]
        end
    end

    if s[1] == 1
        i,j = res[1]
        h = AlgebraElement{C}(lie(s[2]; C), eqs.one)
        eq = (sum([c*kappa[m[1][2],j] for (c, m) in normalForm(eqs.d, comm(h, mod(i; C)))]; init=AlgebraElement{C}()) # κ([h⋅v_i,v_j])
            + sum([c*kappa[i,m[1][2]] for (c, m) in normalForm(eqs.d, comm(h, mod(j; C)))]; init=AlgebraElement{C}()) # κ([v_i,h⋅v_j])
            - normalForm(eqs.d, comm(h, kappa[i,j])))                                                                 # h⋅κ([v_i,v_j])
    else
        i,j,k = res[1]
        eq = (kappa[i,j]*mod(k; C) - mod(i; C)*kappa[j,k]
            + kappa[j,k]*mod(i; C) - mod(j; C)*kappa[k,i]
            + kappa[k,i]*mod(j; C) - mod(k; C)*kappa[i,j]
            - kappa[k,j]*mod(i; C) + mod(k; C)*kappa[j,i]
            - kappa[j,i]*mod(k; C) + mod(j; C)*kappa[i,k]
            - kappa[i,k]*mod(j; C) + mod(i; C)*kappa[k,j])
    end

    return (normalForm(eqs.d, eq), s)
end

function Base.length(eqs::PBWDeformEqs{C}) where C
    nL = eqs.d.extraData.sp.nL
    nV = eqs.d.extraData.sp.nV
    return binomial(nV,2)*nL + binomial(nV,3)
end


function PBWDeformEqs1(d::QuadraticAlgebra{C, SmashProductDeformLie{C}}, one=C(1)::C) :: Vector{AlgebraElement{C}} where C
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras.
    # DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011
    nL = d.extraData.sp.nL
    nV = d.extraData.sp.nV
    kappa = d.extraData.kappa

    ## (a) κ is H-invariant
    eqs = [(sum([c*kappa[m[1][2],j] for (c, m) in normalForm(d, comm(h, mod(i; C)))]; init=AlgebraElement{C}()) # κ([h⋅v_i,v_j])
          + sum([c*kappa[i,m[1][2]] for (c, m) in normalForm(d, comm(h, mod(j; C)))]; init=AlgebraElement{C}()) # κ([v_i,h⋅v_j])
          - normalForm(d, comm(h, kappa[i,j])))                                                                 # h⋅κ([v_i,v_j])
        for h in map(b -> AlgebraElement{C}(b, one), lie(1:nL; C)) for i in 1:nV for j in i+1:nV]
    # m[1][2] denotes the index of the only basis element in the monomial m

    ## (b) trivial

    ## (c) 0 = κ ⊗ id - id ⊗ κ on (I ⊗ V) ∩ (V ⊗ I)
    # (I ⊗ V) ∩ (V ⊗ I) has basis v_iv_jv_k + v_jv_kv_i + v_kv_iv_j - v_kv_jv_i - v_jv_iv_k - v_iv_kv_j for i<j<k
    append!(eqs, [(kappa[i,j]*mod(k; C) - mod(i; C)*kappa[j,k]
            + kappa[j,k]*mod(i; C) - mod(j; C)*kappa[k,i]
            + kappa[k,i]*mod(j; C) - mod(k; C)*kappa[i,j]
            - kappa[k,j]*mod(i; C) + mod(k; C)*kappa[j,i]
            - kappa[j,i]*mod(k; C) + mod(j; C)*kappa[i,k]
            - kappa[i,k]*mod(j; C) + mod(i; C)*kappa[k,j]
        ) for i in 1:nV for j in i+1:nV for k in j+1:nV])

    ## (d) trivial

    return map(eq -> normalForm(d, eq), eqs)
end

function isPBWDeform(d::QuadraticAlgebra{C, SmashProductDeformLie{C}}, one=C(1)::C) :: Bool where C
    return all(iszero, PBWDeformEqs{C}(d, one))
end


function paramDeformNumberVars(nL::Int64, nV::Int64, maxdeg::Int64) :: Tuple{Int64, Int64, Int64}
    nKappaEntries = div(nV*(nV-1), 2)
    nEntryCoeffs = sum(multichoose(nL, k) for k in 0:maxdeg)
    
    return nKappaEntries * nEntryCoeffs, nKappaEntries, nEntryCoeffs
end

function paramDeformVars(nL::Int64, nV::Int64, maxdeg::Int64) :: Vector{String}
    # format: "c_{i,j,deg,[inds]}"
    return ["c_{$i,$j,$deg,$(isempty(inds) ? "[]" : inds)}" for i in 1:nV for j in i+1:nV for deg in 0:maxdeg for inds=multicombinations(1:nL, deg)]
end

function sortVars(vars::Vector{T}, nL, nV, maxdeg) :: Matrix{Vector{Vector{T}}} where T
    _, nKappaEntries, nEntryCoeffs = paramDeformNumberVars(nL, nV, maxdeg)
    m = fill(Vector{T}[], nV, nV)
    k = 0
    for i in 1:nV, j in i+1:nV
        offset = 0
        m[i,j] = fill(T[], maxdeg+1)
        for d in 0:maxdeg
            curr = multichoose(nL, d)
            m[i,j][d+1] = vars[k*nEntryCoeffs+1+offset : k*nEntryCoeffs+offset+curr]
            offset += curr
        end
        k += 1
    end
    return m
end



function varietyOfPBWDeformsLinear(sp::QuadraticAlgebra{Rational{Int64}, SmashProductLie}, maxdeg::Int64) :: SparseArrays.SparseMatrixCSC{fmpq, Int64}
    nL = sp.extraData.nL
    nV = sp.extraData.nV
  
    log("Constructing MPolyRing...")
    R, vars = PolynomialRing(QQ, paramDeformVars(nL, nV, maxdeg))
    numVars = length(vars)
    varLookup = Dict(vars[i] => i for i in 1:numVars)

    varMatrix = sortVars(vars, nL, nV, maxdeg)

    log("Constructing kappa...")
    kappa = fill(AlgebraElement{MPolyElem}(0), nV, nV)
    for i in 1:nV, j in i+1:nV, d in 0:maxdeg, (k, ind) in enumerate(multicombinations(1:nL, d))
        kappa[i,j] += varMatrix[i,j][d+1][k]*lie(ind; C=MPolyElem)
        kappa[j,i] -= varMatrix[i,j][d+1][k]*lie(ind; C=MPolyElem)
    end

    log("Changing SmashProductLie coeffcient type...")
    newBasis = [changeC(MPolyElem, b) for b in sp.basis]
    newRelTable = Dict([(changeC(MPolyElem, b1), changeC(MPolyElem, b2)) => 
        AlgebraElement{MPolyElem}(map(x -> (R(x[1]), changeC(MPolyElem, x[2])), unpack(a))) 
        for ((b1, b2), a) in pairs(sp.relTable)])
    newSp = QuadraticAlgebra{MPolyElem, SmashProductLie}(newBasis, newRelTable, sp.extraData)

    log("Constructing deformation...")
    deform = smashProductDeformLie(newSp, kappa, R(1))

    log("Generating equation iterator...")
    iter = Iterators.map(a -> poly2vecLinear(a, varLookup, numVars),
        Iterators.flatten(
            Iterators.map(coefficientComparison,
                PBWDeformEqs{MPolyElem}(deform, R(1))
            )
        )
    )

    # group sparse vectors by index of first non-zero entry
    log("Collecting rows...")
    lgs = [Vector{SparseVector{fmpq, Int64}}() for _ in 1:numVars]
    for v in iter
        normalizeAndStore!(lgs, v)
    end

    # create row-echelon form
    log("Computing row-echelon form...")
    row_echelon!(lgs)

    # reduce row-echelon form
    log("Computing reduced row-echelon form...")
    reduced_row_echelon!(lgs)

    return lgs2mat(lgs, numVars)
end

@inline function normalizeAndStore!(lgs::Vector{Vector{SparseVector{T, Int64}}}, v::SparseVector{T, Int64}) where {T <: AbstractAlgebra.RingElement}
    nzIndices, nzValues = findnz(v)
    push!(lgs[nzIndices[1]], inv(nzValues[1]) .* v)
end

@inline function poly2vecLinear(a::fmpq_mpoly, varLookup::Dict{fmpq_mpoly, Int64}, numVars::Int64) :: SparseVector{fmpq, Int64}
    @assert total_degree(a) == 1

    return sparsevec(
        Dict(varLookup[monomial(a,i)] => coeff(a,i) for i in 1:length(a)),
        numVars
    )
end

function row_echelon!(lgs::Vector{Vector{SparseVector{T, Int64}}}) :: Vector{Vector{SparseVector{T, Int64}}} where {T <: AbstractAlgebra.RingElement}
    for i in 1:length(lgs)
        unique!(lgs[i])
        if length(lgs[i]) <= 1
            continue
        end

        for j in 2:length(lgs[i])
            lgs[i][j] -= lgs[i][1]
            if !iszero(lgs[i][j])
                normalizeAndStore!(lgs, lgs[i][j])
            end
        end
        deleteat!(lgs[i], 2:length(lgs[i]))
    end
    return lgs
end

function reduced_row_echelon!(lgs::Vector{Vector{SparseVector{T, Int64}}}) :: Vector{Vector{SparseVector{T, Int64}}} where {T <: AbstractAlgebra.RingElement}
    for i in length(lgs):-1:1
        if isempty(lgs[i])
            continue
        end
        nzIndices, nzValues = findnz(lgs[i][1])
        for (ind,j) in enumerate(nzIndices[2:end])
            if !isempty(lgs[j])
                lgs[i][1] -= nzValues[ind+1] .* lgs[j][1]
            end
        end
    end
    return lgs
end

function lgs2mat(lgs::Vector{Vector{SparseVector{T, Int64}}}, n::Int64) :: SparseArrays.SparseMatrixCSC{T, Int64}  where {T <: AbstractAlgebra.RingElement}
    mat = spzeros(fmpq, n, n)
    for i in 1:n
        if !isempty(lgs[i])
            mat[i,:] = lgs[i][1]
        end
    end
    return mat
end

function coefficientComparison(eq::AlgebraElement{C}) :: Vector{C} where C
    result = C[]
    for summand in unpack(eq)
        (c, m) = summand
        push!(result, c)
    end
    return result
end

function simplifyGen(gen::MPolyElem) :: MPolyElem
    return ((leading_coefficient(gen) < 0 ? -1 : 1)//content(gen))*gen
end


function varietyOfPBWDeforms1(sp::QuadraticAlgebra{Rational{Int64}, SmashProductLie}, maxdeg::Int64) :: SparseArrays.SparseMatrixCSC{fmpq, Int64}
    nL = sp.extraData.nL
    nV = sp.extraData.nV

    log("Constructing MPolyRing...")
    R, vars = PolynomialRing(QQ, paramDeformVars(nL, nV, maxdeg))
    numVars = length(vars)
    varLookup = Dict(vars[i] => i for i in 1:numVars)
    varMatrix = sortVars(vars, nL, nV, maxdeg)

    log("Constructing kappa...")
    kappa = fill(AlgebraElement{MPolyElem}(0), nV, nV)
    for i in 1:nV, j in i+1:nV, d in 0:maxdeg, (k, ind) in enumerate(multicombinations(1:nL, d))
        kappa[i,j] += varMatrix[i,j][d+1][k]*lie(ind; C=MPolyElem)
        kappa[j,i] -= varMatrix[i,j][d+1][k]*lie(ind; C=MPolyElem)
    end

    log("Changing SmashProductLie coeffcient type...")
    newBasis = [changeC(MPolyElem, b) for b in sp.basis]
    newRelTable = Dict([(changeC(MPolyElem, b1), changeC(MPolyElem, b2)) =>
        AlgebraElement{MPolyElem}(map(x -> (R(x[1]), changeC(MPolyElem, x[2])), unpack(a)))
        for ((b1, b2), a) in pairs(sp.relTable)])
    newSp = QuadraticAlgebra{MPolyElem, SmashProductLie}(newBasis, newRelTable, sp.extraData)

    log("Constructing deformation...")
    deform = smashProductDeformLie(newSp, kappa, R(1))

    log("Generating equations...")
    iter = map(a -> poly2vecLinear(a, varLookup, numVars), coefficientComparison1(PBWDeformEqs1(deform, R(1))))

    # group sparse vectors by index of first non-zero entry
    log("Collecting rows...")
    lgs = [Vector{SparseVector{fmpq, Int64}}() for _ in 1:numVars]
    for v in iter
        normalizeAndStore!(lgs, v)
    end

    # create row-echelon form
    log("Computing row-echelon form...")
    row_echelon!(lgs)

    # reduce row-echelon form
    log("Computing reduced row-echelon form...")
    reduced_row_echelon!(lgs)

    return lgs2mat(lgs, numVars)
end

function coefficientComparison1(eqs::Vector{AlgebraElement{C}}) :: Vector{C} where C
    result = C[]
    for eq in eqs
        # coeffcient comparison
        for summand in unpack(eq)
            (c, m) = summand
            push!(result, c)
        end
    end
    return result
end

function simplifyGens1(gens::Vector{C}) :: Vector{C} where C
    gens = [(1//content(gen))*gen for gen in gens]
    gens = [(leading_coefficient(gen) < 0 ? -gen : gen) for gen in gens]
    return unique(gens)
end
