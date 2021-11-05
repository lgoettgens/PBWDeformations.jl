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

    free_alg = free_algebra(coeff_ring, [symbL; symbV])
    free_baseL = [gen(free_alg, i) for i in 1:dimL]
    free_baseV = [gen(free_alg, dimL+i) for i in 1:dimV]

    rels = Dict{Tuple{Int,Int}, FreeAlgebraElem{C}}()

    for i in 1:dimL, j in 1:dimL
        rels[(i, j)] = free_baseL[j] * free_baseL[i] + sum(c * free_baseL[k] for (c, k) in struct_const_L[i,j]; init=0)
    end

    for i in 1:dimL, j in 1:dimV
        rels[(i, dimL+j)] = free_baseV[j] * free_baseL[i] + sum(c * free_baseV[k] for (c, k) in struct_const_V[i,j]; init=0)
        rels[(dimL+j, i)] = -rels[(i, dimL+j)]
    end

    alg = quadratic_quo_algebra(free_alg, rels)
    baseL = [gen(alg, i) for i in 1:dimL]
    baseV = [gen(alg, dimL+i) for i in 1:dimV]

    return SmashProductLie{C}(dimL, dimV, baseL, baseV, coeff_ring, alg), (baseL, baseV)
end

function smash_product_lie(coeff_ring :: Ring, symbL :: Vector{String}, symbV :: Vector{String}, struct_const_L :: Matrix{Vector{Tuple{Int, Int}}}, struct_const_V :: Matrix{Vector{Tuple{Int, Int}}})
    return smash_product_lie(coeff_ring, map(Symbol, symbL), map(Symbol, symbV), struct_const_L, struct_const_V)
end

function smash_product_lie(coeff_ring :: Ring, dynkin :: Char, n :: Int, lambda :: Vector{Int})
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

    symbL = ["x_$i" for i in 1:dimL]
    symbV = ["v_$i" for i in 1:dimV]

    return smash_product_lie(coeff_ring, symbL, symbV, struct_const_L, struct_const_V)
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
