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


function PBWDeformEqs(d::QuadraticAlgebra{C, SmashProductDeformLie{C}}, one=C(1)::C) :: Vector{AlgebraElement{C}} where C
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras.
    # DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011
    nL = d.extraData.sp.nL
    nV = d.extraData.sp.nV
    kappa = d.extraData.kappa

    ## (a) κ is H-invariant
    eqs = [(sum([c*kappa[m[1][2],j] for (c, m) in normalForm(d, comm(h, mod(i; C)))]; init=AlgebraElement{C}()) # κ([h⋅v_i,v_j])
          + sum([c*kappa[i,m[1][2]] for (c, m) in normalForm(d, comm(h, mod(j; C)))]; init=AlgebraElement{C}()) # κ([v_i,h⋅v_j])
          - normalForm(d, comm(h, kappa[i,j])))                                                                 # h⋅κ([v_i,v_j])
        for i in 1:nV for j in i+1:nV for h in map(b -> AlgebraElement{C}(b, one), lie(1:nL; C))]
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
    return all(iszero, PBWDeformEqs(d, one))
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

function varietyOfPBWDeforms(sp::QuadraticAlgebra{Rational{Int64}, SmashProductLie}, maxdeg::Int64) :: Tuple{Vector{MPolyElem}, MPolyRing}
    nL = sp.extraData.nL
    nV = sp.extraData.nV
  

    R, vars = PolynomialRing(QQ, paramDeformVars(nL, nV, maxdeg))

    varMatrix = sortVars(vars, nL, nV, maxdeg)

    kappa = fill(AlgebraElement{MPolyElem}(0), nV, nV)
    for i in 1:nV, j in i+1:nV, d in 0:maxdeg, (k, ind) in enumerate(multicombinations(1:nL, d))
        kappa[i,j] += varMatrix[i,j][d+1][k]*lie(ind; C=MPolyElem)
        kappa[j,i] -= varMatrix[i,j][d+1][k]*lie(ind; C=MPolyElem)
    end

    newBasis = [changeC(MPolyElem, b) for b in sp.basis]
    newRelTable = Dict([(changeC(MPolyElem, b1), changeC(MPolyElem, b2)) => 
        AlgebraElement{MPolyElem}(map(x -> (R(x[1]), changeC(MPolyElem, x[2])), unpack(a))) 
        for ((b1, b2), a) in pairs(sp.relTable)])
    newSp = QuadraticAlgebra{MPolyElem, SmashProductLie}(newBasis, newRelTable, sp.extraData)

    deform = smashProductDeformLie(newSp, kappa, R(1))

    return simplifyGens(coefficientComparison(PBWDeformEqs(deform, R(1)))), R
end

function coefficientComparison(eqs::Vector{AlgebraElement{C}}) :: Vector{C} where C
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

function simplifyGens(gens::Vector{C}) :: Vector{C} where C
    gens = [(1//content(gen))*gen for gen in gens]
    gens = [(leading_coefficient(gen) < 0 ? -gen : gen) for gen in gens]
    return unique(gens)
end
