mutable struct SmashProductLie{C <: RingElement}
    dimL :: Int64
    dimV :: Int64
    baseL :: Vector{QuadraticQuoAlgebraElem{C}}
    baseV :: Vector{QuadraticQuoAlgebraElem{C}}
    coeff_ring :: Ring
    alg :: QuadraticQuoAlgebra{C}
    # dynkin :: Char
    # n :: Int64
    # lambda :: Vector{Int64}
    # matrixRepL :: Vector{Matrix{Int64}}
end


function smash_product_lie(coeff_ring :: Ring, symbL :: Vector{Symbol}, symbV :: Vector{Symbol}, struct_const_L :: Matrix{Vector{Tuple{Int, Int}}}, struct_const_V :: Matrix{Vector{Tuple{Int, Int}}})
    C = elem_type(coeff_ring)
    
    dimL = length(symbL)
    dimV = length(symbV)

    free_alg, _ = free_algebra(coeff_ring, [symbL; symbV])
    free_baseL = [gen(free_alg, i) for i in 1:dimL]
    free_baseV = [gen(free_alg, dimL+i) for i in 1:dimV]

    rels = Dict{Tuple{Int,Int}, FreeAlgebraElem{C}}()

    for i in 1:dimL, j in 1:dimL
        rels[(i, j)] = free_baseL[j] * free_baseL[i] + sum(c * free_baseL[k] for (c, k) in struct_const_L[i,j]; init=zero(free_alg))
    end

    for i in 1:dimL, j in 1:dimV
        rels[(i, dimL+j)] = free_baseV[j] * free_baseL[i] + sum(c * free_baseV[k] for (c, k) in struct_const_V[i,j]; init=zero(free_alg))
        rels[(dimL+j, i)] = free_baseL[i] * free_baseV[j] - sum(c * free_baseV[k] for (c, k) in struct_const_V[i,j]; init=zero(free_alg))
    end

    alg, _ = quadratic_quo_algebra(free_alg, rels)
    baseL = [gen(alg, i) for i in 1:dimL]
    baseV = [gen(alg, dimL+i) for i in 1:dimV]

    return SmashProductLie{C}(dimL, dimV, baseL, baseV, coeff_ring, alg), (baseL, baseV)
end

function smash_product_lie(coeff_ring :: Ring, symbL :: Vector{String}, symbV :: Vector{String}, struct_const_L :: Matrix{Vector{Tuple{Int, Int}}}, struct_const_V :: Matrix{Vector{Tuple{Int, Int}}})
    return smash_product_lie(coeff_ring, map(Symbol, symbL), map(Symbol, symbV), struct_const_L, struct_const_V)
end

function smash_product_lie(coeff_ring :: Ring, dynkin :: Char, n :: Int, lambda :: Vector{Int})
    dimL, dimV, struct_const_L, struct_const_V = smash_product_struct_const_from_gap(dynkin, n, lambda)

    symbL = ["x_$i" for i in 1:dimL]
    symbV = ["v_$i" for i in 1:dimV]

    return smash_product_lie(coeff_ring, symbL, symbV, struct_const_L, struct_const_V)
end

function smash_product_struct_const_from_gap(dynkin :: Char, n :: Int, lambda :: Vector{Int})
    @assert n == length(lambda)
    sanitize_lie_input(dynkin, n)

    GAPG = GAP.Globals

    L = GAPG.SimpleLieAlgebra(GAP.julia_to_gap(string(dynkin)), n, GAPG.Rationals)
    dimL = GAPG.Dimension(L)
    basisL = GAPG.BasisVectors(GAPG.Basis(L))
    comm_table_L = GAP.gap_to_julia(GAPG.StructureConstantsTable(GAPG.Basis(L)))[1:dimL]

    struct_const_L = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimL)
    for i in 1:dimL, j in 1:dimL
        struct_const_L[i,j] = [(c, k) for (k, c) in zip(comm_table_L[i][j]...)]
    end

    V = GAPG.HighestWeightModule(L, GAP.julia_to_gap(lambda))
    dimV = GAPG.Dimension(V)
    basisV = GAPG.BasisVectors(GAPG.Basis(V))
    
    struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
    for i in 1:dimL, j in 1:dimV
        struct_const_V[i,j] = [(c, k) for (k, c) in enumerate(GAP.gap_to_julia(GAPG.Coefficients(GAPG.Basis(V), basisL[i]^basisV[j]))) if !iszero(c)]
    end

    return dimL, dimV, struct_const_L, struct_const_V
end

function smash_product_struct_const_so(n :: Int, lambda :: Vector{Int}) # for odd n this is a B type highest weight, for even n D type
    @assert div(n,2) == length(lambda)
    lenghtlambda = div(n,2)

    ur_triag(M) = vcat([M[i, i+1:end] for i in 1:size(M,1)]...)
    std_basis(i,n) = [i == j ? 1 : 0 for j in 1:n]

    dimL = div(n*(n-1), 2)
    basisL = [(b = zeros(Int,n,n); b[i,j] = 1; b[j,i] = -1; b) for i in 1:n for j in i+1:n]
    symbL = ["x_$(i)_$(j)" for i in 1:n for j in i+1:n]

    struct_const_L = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimL)
    for i in 1:dimL, j in 1:dimL
        struct_const_L[i,j] = [(c, k) for (k, c) in enumerate(ur_triag(basisL[i]*basisL[j] - basisL[j]*basisL[i])) if !iszero(c)]
    end

    if in(lambda, [std_basis(i,lenghtlambda) for i in 1:lenghtlambda])
        lambda_std_i = findfirst(==(lambda), [std_basis(i,lenghtlambda) for i in 1:lenghtlambda])

        if (n % 2 == 1 && lambda_std_i <= lenghtlambda-1) || (n % 2 == 0 && lambda_std_i <= lenghtlambda-2) # exterior product of defining representation
            if lambda_std_i == 1
                dimV = n
                symbV = ["v_$(i)" for i in 1:dimV]
                struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
                for i in 1:dimL, j in 1:dimV
                    struct_const_V[i,j] = [(c, k) for (k, c) in enumerate(basisL[i] * std_basis(j,n)) if !iszero(c)]
                end
            elseif lambda_std_i == 2
                dimV = binomial(n,2)
                symbV = ["v_$(i)_$(j)" for i in 1:n for j in i+1:n]
                struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
                for i in 1:dimL, j1 in 1:dimV, j2 in j1+1:n
                    baseIndex(j1,j2) = binomial(n,2) - binomial(n-j1+1,2) + j2 - j1
                    struct_const_V[i,baseIndex(j1,j2)] = [
                        collect(k1 < j2 ? (c, baseIndex(k1,j2)) : (-c, baseIndex(j2,k1)) for (k1, c) in enumerate(basisL[i] * std_basis(j1,n)) if !iszero(c) && k1 != j2);
                        collect(j1 < k2 ? (c, baseIndex(j1,k2)) : (-c, baseIndex(k2,j1)) for (k2, c) in enumerate(basisL[i] * std_basis(j2,n)) if !iszero(c) && j1 != k2)
                    ]
                end
            else
                error("Not implemented yet.")
            end
        else
            error("Spin representations are not implemented yet.")
        end
    else
        error("Non-fundamental representations are not implemented yet.")
    end

    return dimL, dimV, struct_const_L, struct_const_V, symbL, symbV
end


ngens(sp::SmashProductLie) = sp.dimL, sp.dimV

function gens(sp::SmashProductLie{C}) where C <: RingElement
    return [gen(sp.alg, i) for i in 1:sp.dimL], [gen(sp.alg, i+sp.dimL) for i in 1:sp.dimV]
end


function show(io::IO, sp::SmashProductLie)
    local max_gens = 4 # largest number of generators to print
    print(io, "Lie Algebra Smash Product with basis ")
    for i = 1:min(sp.dimL - 1, max_gens - 1)
       print(io, string(sp.alg.S[i]), ", ")
    end
    if sp.dimL > max_gens
       print(io, "..., ")
    end
    print(io, string(sp.alg.S[sp.dimL]) * ", ")
    for i = 1:min(sp.dimV - 1, max_gens - 1)
        print(io, string(sp.alg.S[sp.dimL+i]), ", ")
     end
     if sp.dimV > max_gens
        print(io, "..., ")
     end
     print(io, string(sp.alg.S[sp.dimL+sp.dimV]))
    print(io, " over ")
    print(IOContext(io, :compact => true), sp.coeff_ring)
end


function change_base_ring(R::Ring, sp::SmashProductLie{C}) where C <: RingElement
    alg = change_base_ring(R, sp.alg)
    baseL = [gen(alg, i) for i in 1:sp.dimL]
    baseV = [gen(alg, sp.dimL+i) for i in 1:sp.dimV]

    return SmashProductLie{elem_type(R)}(sp.dimL, sp.dimV, baseL, baseV, R, alg)
end
