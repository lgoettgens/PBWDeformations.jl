struct SmashProductDeformLie{C <: ScalarTypes}
    sp :: SmashProductLie
    symmetric :: Bool
    kappa :: Matrix{AlgebraElement{C}}
end

function Base.:(==)(spd1::SmashProductDeformLie{C}, spd2::SmashProductDeformLie{C}) :: Bool where C <: ScalarTypes
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


function smash_product_deform_lie(sp::QuadraticAlgebra{C, SmashProductLie}, kappa::Matrix{AlgebraElement{C}}, one=C(1)::C) :: QuadraticAlgebra{C, SmashProductDeformLie{C}} where C <: ScalarTypes
    nV = sp.extraData.nV
    @assert size(kappa) == (nV, nV) "size of kappa matches module dimension"

    # basis of smash product consists of basis of module and basis of Hopf algebra
    hopfBasis = filter(!ismod, sp.basis)
    @assert all(e -> issubset(basis_elements(e), hopfBasis), kappa) "kappa only takes values in Hopf algebra"

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


function smash_product_symm_deform_lie(sp::QuadraticAlgebra{C, SmashProductLie}) :: QuadraticAlgebra{C, SmashProductDeformLie{C}} where C <: ScalarTypes
    relTable = sp.relTable

    for i in 1:sp.extraData.nV, j in 1:i-1
        relTable[(mod(i; C), mod(j; C))] = AlgebraElement{C}(mod(j, i; C))
    end

    extraData = SmashProductDeformLie{C}(sp.extraData, true, fill(AlgebraElement{C}(), sp.extraData.nV, sp.extraData.nV))
    return QuadraticAlgebra{C, SmashProductDeformLie{C}}(sp.basis, relTable, extraData)
end

function smash_product_symm_deform_lie(dynkin::Char, n::Int64, lambda::Vector{Int64}; C::Type{<:ScalarTypes} = DefaultScalarType) :: QuadraticAlgebra{C, SmashProductDeformLie{C}}
    @assert n == length(lambda)
    sanitize_lie_input(dynkin, n)

    return smash_product_symm_deform_lie(smash_product_lie(dynkin, n, lambda; C))
end

struct PBWDeformEqs{C <: ScalarTypes}
    d :: QuadraticAlgebra{C, SmashProductDeformLie{C}}
    one :: C

    PBWDeformEqs{C}(d::QuadraticAlgebra{C, SmashProductDeformLie{C}}, one = C(1)) where C <: ScalarTypes = new{C}(d, one)
end

function Base.iterate(eqs::PBWDeformEqs{C}, s::Tuple{Int64, Int64, Union{Nothing, Vector{Int64}}} = (0, 0, nothing)) where C <: ScalarTypes
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras.
    # DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011

    nL = eqs.d.extraData.sp.nL
    nV = eqs.d.extraData.sp.nV
    kappa = eqs.d.extraData.kappa

    # The following rather complicated code computes the parameters for the current equation and the next state, using the given state of the previous iteration
    #s = (phase = 1, h, c_state) or (phase = 3, counter, c_state)

    if nV < 2 # for these cases no equations exist
        return nothing
    end

    if s[1] == 0
        s = (1, 1, nothing)
    end

    if s[1] == 1
        comb = Combinatorics.Combinations(nV,2)
        res = s[3] === nothing ? iterate(comb) : iterate(comb, s[3])
        if res === nothing
            @debug "Equation generation, first phase $(floor(Int, 100*s[2] / nL))%"
            s = (s[1], s[2]+1, s[3])
            if s[2] > nL
                s = (3, 0, nothing)
            else
                comb = Combinatorics.Combinations(nV,2)
                res = iterate(comb)
                s = (s[1], s[2], res[2])
            end
        else
            s = (s[1], s[2], res[2])
        end
    end
    if s[1] == 3
        comb = Combinatorics.Combinations(nV,3)
        res = s[3] === nothing ? iterate(comb) : iterate(comb, s[3])
        if res === nothing
            return nothing
        else
            s = (s[1], s[2], res[2])
        end
        s = (s[1], s[2]+1, s[3])
        if (s[2] % 100 == 0)
            @debug "Equation generation, second phase $(floor(Int, 100*s[2] / binomial(nV, 3)))%"
        end
    end

    ## (a) κ is H-invariant
    if s[1] == 1
        i,j = res[1]
        h = AlgebraElement{C}(lie(s[2]; C), eqs.one)
        eq = (sum([c*kappa[m[1][2],j] for (c, m) in normal_form(eqs.d, comm(h, mod(i; C)))]; init=AlgebraElement{C}()) # κ([h⋅v_i,v_j])
            + sum([c*kappa[i,m[1][2]] for (c, m) in normal_form(eqs.d, comm(h, mod(j; C)))]; init=AlgebraElement{C}()) # κ([v_i,h⋅v_j])
            - normal_form(eqs.d, comm(h, kappa[i,j])))                                                                 # h⋅κ([v_i,v_j])
        # m[1][2] denotes the index of the only basis element in the monomial m

    ## (b) trivial

    ## (c) 0 = κ ⊗ id - id ⊗ κ on (I ⊗ V) ∩ (V ⊗ I)
    # (I ⊗ V) ∩ (V ⊗ I) has basis v_iv_jv_k + v_jv_kv_i + v_kv_iv_j - v_kv_jv_i - v_jv_iv_k - v_iv_kv_j for i<j<k
    elseif s[1] == 3
        i,j,k = res[1]
        eq = (kappa[i,j]*mod(k; C) - mod(i; C)*kappa[j,k]
            + kappa[j,k]*mod(i; C) - mod(j; C)*kappa[k,i]
            + kappa[k,i]*mod(j; C) - mod(k; C)*kappa[i,j]
            - kappa[k,j]*mod(i; C) + mod(k; C)*kappa[j,i]
            - kappa[j,i]*mod(k; C) + mod(j; C)*kappa[i,k]
            - kappa[i,k]*mod(j; C) + mod(i; C)*kappa[k,j])

    ## (d) trivial
    else
        @error "This should nod be reached."
    end

    return (normal_form(eqs.d, eq), s)
end

function Base.length(eqs::PBWDeformEqs{C}) where C <: ScalarTypes
    nL = eqs.d.extraData.sp.nL
    nV = eqs.d.extraData.sp.nV
    return binomial(nV,2)*nL + binomial(nV,3)
end


function pbwdeform_eqs_noiter(d::QuadraticAlgebra{C, SmashProductDeformLie{C}}, one = C(1)::C) :: Vector{AlgebraElement{C}} where C <: ScalarTypes
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras.
    # DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011
    nL = d.extraData.sp.nL
    nV = d.extraData.sp.nV
    kappa = d.extraData.kappa

    @debug "Equation generation, first phase"

    ## (a) κ is H-invariant
    eqs = [(sum([c*kappa[m[1][2],j] for (c, m) in normal_form(d, comm(h, mod(i; C)))]; init=AlgebraElement{C}()) # κ([h⋅v_i,v_j])
          + sum([c*kappa[i,m[1][2]] for (c, m) in normal_form(d, comm(h, mod(j; C)))]; init=AlgebraElement{C}()) # κ([v_i,h⋅v_j])
          - normal_form(d, comm(h, kappa[i,j])))                                                                 # h⋅κ([v_i,v_j])
        for h in map(b -> AlgebraElement{C}(b, one), lie(1:nL; C)) for i in 1:nV for j in i+1:nV]
    # m[1][2] denotes the index of the only basis element in the monomial m

    ## (b) trivial

    @debug "Equation generation, second phase"

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

    return map(eq -> normal_form(d, eq), eqs)
end

function ispbwdeform(d::QuadraticAlgebra{C, SmashProductDeformLie{C}}, one = C(1)::C) :: Bool where C <: ScalarTypes
    return all(iszero, PBWDeformEqs{C}(d, one))
end


function param_deform_number_vars(nL::Int64, nV::Int64, maxdeg::Int64) :: Tuple{Int64, Int64, Int64}
    nKappaEntries = div(nV*(nV-1), 2)
    nEntryCoeffs = sum(binomial(nL + k - 1, k) for k in 0:maxdeg)
    
    return nKappaEntries * nEntryCoeffs, nKappaEntries, nEntryCoeffs
end

function param_deform_vars(nL::Int64, nV::Int64, maxdeg::Int64) :: Vector{String}
    # format: "c_{i,j,deg,[inds]}"
    return ["c_{$i,$j,$deg,$(isempty(inds) ? "[]" : inds)}" for i in 1:nV for j in i+1:nV for deg in 0:maxdeg for inds=Combinatorics.with_replacement_combinations(1:nL, deg)]
end

function sort_vars(vars::Vector{T}, nL, nV, maxdeg) :: Matrix{Vector{Vector{T}}} where T
    _, nKappaEntries, nEntryCoeffs = param_deform_number_vars(nL, nV, maxdeg)
    m = fill(Vector{T}[], nV, nV)
    k = 0
    for i in 1:nV, j in i+1:nV
        offset = 0
        m[i,j] = fill(T[], maxdeg+1)
        for d in 0:maxdeg
            curr = binomial(nL + d - 1, d)
            m[i,j][d+1] = vars[k*nEntryCoeffs+1+offset : k*nEntryCoeffs+offset+curr]
            offset += curr
        end
        k += 1
    end
    return m
end


function possible_pbw_deforms(sp::QuadraticAlgebra{DefaultScalarType, SmashProductLie}, maxdeg::Int64;
            use_iterators=true::Bool, special_return::Type{T}=Nothing) where T <: Union{Nothing, SparseMatrixCSC}

    nL = sp.extraData.nL
    nV = sp.extraData.nV

    @info "Constructing MPolyRing..."
    R, vars = PolynomialRing(QQ, param_deform_vars(nL, nV, maxdeg))
    numVars = length(vars)
    varLookup = Dict(vars[i] => i for i in 1:numVars)
    varMatrix = sort_vars(vars, nL, nV, maxdeg)

    @info "Constructing kappa..."
    kappa = fill(AlgebraElement{MPolyElem}(0), nV, nV)
    for i in 1:nV, j in i+1:nV, d in 0:maxdeg, (k, ind) in enumerate(Combinatorics.with_replacement_combinations(1:nL, d))
        kappa[i,j] += varMatrix[i,j][d+1][k]*lie(ind; C=MPolyElem)
        kappa[j,i] -= varMatrix[i,j][d+1][k]*lie(ind; C=MPolyElem)
    end

    @info "Changing SmashProductLie coeffcient type..."
    newBasis = [change_c(MPolyElem, b) for b in sp.basis]

    newRelTable = Dict([(change_c(MPolyElem, b1), change_c(MPolyElem, b2)) =>
        AlgebraElement{MPolyElem}(map(x -> (R(x[1]), change_c(MPolyElem, x[2])), unpack(a)))
        for ((b1, b2), a) in pairs(sp.relTable)])
    newSp = QuadraticAlgebra{MPolyElem, SmashProductLie}(newBasis, newRelTable, sp.extraData)

    @info "Constructing deformation..."
    deform = smash_product_deform_lie(newSp, kappa, R(1))

    if use_iterators
        @info "Generating equation iterator..."
        iter = Iterators.map(a -> poly2vec_linear(a, varLookup, numVars),
            Iterators.flatten(
                Iterators.map(coefficient_comparison,
                    PBWDeformEqs{MPolyElem}(deform, R(1))
                )
            )
        )
    else
        @info "Generating equations..."
        iter = map(a -> poly2vec_linear(a, varLookup, numVars), reduce(vcat, map(coefficient_comparison, pbwdeform_eqs_noiter(deform, R(1)))))
    end

    # group sparse vectors by index of first non-zero entry
    @info "Collecting rows..."
    lgs = [Vector{SparseVector{fmpq, Int64}}() for _ in 1:numVars]
    for v in iter
        normalize_and_store!(lgs, v)
    end

    # create row-echelon form
    @info "Computing row-echelon form..."
    row_echelon!(lgs)

    # reduce row-echelon form
    @info "Computing reduced row-echelon form..."
    reduced_row_echelon!(lgs)

    mat = lgs2mat(lgs, numVars)
    #return mat, R, vars

    if special_return === SparseMatrixCSC
        return mat, vars
    end

    freedom_ind = indices_of_freedom(mat)
    kappa = fill(AlgebraElement{MPolyElem}(0), nV, nV)
    if length(freedom_ind) > 0
        S, free_params = PolynomialRing(QQ, ["t_$i" for i in 1:length(freedom_ind)])
        for i in 1:nV, j in i+1:nV, d in 0:maxdeg, (k, ind) in enumerate(Combinatorics.with_replacement_combinations(1:nL, d))
            var_ind = varLookup[varMatrix[i,j][d+1][k]]
            if iszero(mat[var_ind,var_ind])
                kappa[i,j] += free_params[findfirst(isequal(var_ind), freedom_ind)] * lie(ind; C=MPolyElem)
                kappa[j,i] -= free_params[findfirst(isequal(var_ind), freedom_ind)] * lie(ind; C=MPolyElem)
            else
                for col in var_ind+1:numVars
                    if !iszero(mat[var_ind,col])
                        kappa[i,j] += -mat[var_ind,col]*free_params[findfirst(isequal(col), freedom_ind)] * lie(ind; C=MPolyElem)
                        kappa[j,i] -= -mat[var_ind,col]*free_params[findfirst(isequal(col), freedom_ind)] * lie(ind; C=MPolyElem)
                    end
                end
            end
        end
        return kappa, free_params
    else
        return kappa, fmpq_mpoly[]
    end
end

function indices_of_freedom(mat::SparseArrays.SparseMatrixCSC{fmpq, Int64}) :: Vector{Int64}
    size(mat)[1] == size(mat)[2] || throw(ArgumentError("Matrix needs to be square."))
    return filter(i -> mat[i,i] == 0, 1:size(mat)[1])
end

@inline function normalize_and_store!(lgs::Vector{Vector{SparseVector{C, Int64}}}, v::SparseVector{C, Int64}) where C <: ScalarTypes
    nzIndices, nzValues = findnz(v)
    push!(lgs[nzIndices[1]], inv(nzValues[1]) .* v)
end

@inline function poly2vec_linear(a::fmpq_mpoly, varLookup::Dict{fmpq_mpoly, Int64}, numVars::Int64) :: SparseVector{fmpq, Int64}
    @assert total_degree(a) == 1

    return sparsevec(
        Dict(varLookup[monomial(a,i)] => coeff(a,i) for i in 1:length(a)),
        numVars
    )
end

function row_echelon!(lgs::Vector{Vector{SparseVector{C, Int64}}}) :: Vector{Vector{SparseVector{C, Int64}}} where C <: ScalarTypes
    for i in 1:length(lgs)
        if (i % 10 == 0)
            @debug "Row echelon, $i/$(length(lgs)), $(floor(Int, 100*i / length(lgs)))%"
        end
        unique!(lgs[i])
        if length(lgs[i]) <= 1
            continue
        end

        for j in 2:length(lgs[i])
            lgs[i][j] -= lgs[i][1]
            if !iszero(lgs[i][j])
                normalize_and_store!(lgs, lgs[i][j])
            end
        end
        deleteat!(lgs[i], 2:length(lgs[i]))
    end
    return lgs
end

function reduced_row_echelon!(lgs::Vector{Vector{SparseVector{C, Int64}}}) :: Vector{Vector{SparseVector{C, Int64}}} where C <: ScalarTypes
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

function lgs2mat(lgs::Vector{Vector{SparseVector{C, Int64}}}, n::Int64) :: SparseArrays.SparseMatrixCSC{C, Int64}  where C <: ScalarTypes
    mat = spzeros(fmpq, n, n)
    for i in 1:n
        if !isempty(lgs[i])
            mat[i,:] = lgs[i][1]
        end
    end
    return mat
end

function coefficient_comparison(eq::AlgebraElement{C}) :: Vector{C} where C <: ScalarTypes
    result = C[]
    for summand in unpack(eq)
        (c, m) = summand
        push!(result, c)
    end
    return result
end
