###############################################################################
#
#       liealgebra_so, i.e. types B and D
#
###############################################################################

"""
    liealgebra_so_basis(n::Int)

Returns a basis of the orthogonal Lie algebra of dimension `n` ``\\mathfrak{so}_n``.

It consists of the matrices ``X_{ij} = E_{ij} - E_{ji}`` with ``i < j`` in row-major order.
"""
function liealgebra_so_basis(n::Int, R::Ring)
    return basis(special_orthogonal_liealgebra(R, n))
end

function liealgebra_so_symbols(n::Int, R::Ring)
    L = special_orthogonal_liealgebra(R, n)
    return symbols(L)
end

function liealgebra_so_struct_const(n::Int, R::Ring) # so_n
    L = special_orthogonal_liealgebra(R, n)
    dimL = ngens(L)

    struct_const_L = Matrix{Vector{Tuple{elem_type(R), Int}}}(undef, dimL, dimL)
    for (i, xi) in enumerate(gens(L)), (j, xj) in enumerate(gens(L))
        struct_const_L[i, j] = [(c, k) for (k, c) in enumerate(_matrix(bracket(xi, xj))) if !iszero(c)]
    end

    return struct_const_L # dimL, struct_const_L
end


function liealgebra_so_symmpowers_standard_module_symbols(n::Int, e::Int, R::Ring)
    L = special_orthogonal_liealgebra(R, n)
    V = symmetric_power(standard_module(L), e)
    return symbols(V)
end

function liealgebra_so_symmpowers_standard_module_struct_const(n::Int, e::Int, R::Ring) # so_n, e-th symm power of the standard rep
    L = special_orthogonal_liealgebra(R, n)
    dimL = ngens(L)

    V = symmetric_power(standard_module(L), e)
    dimV = ngens(V)
    struct_const_V = Matrix{Vector{Tuple{elem_type(R), Int}}}(undef, dimL, dimV)

    for (i, xi) in enumerate(gens(L)), (j, vj) in enumerate(gens(V))
        struct_const_V[i, j] = [(c, k) for (k, c) in enumerate(_matrix(xi * vj)) if !iszero(c)]
    end

    return struct_const_V
end


function liealgebra_so_extpowers_standard_module_symbols(n::Int, e::Int, R::Ring)
    L = special_orthogonal_liealgebra(R, n)
    V = exterior_power(standard_module(L), e)
    return symbols(V)
end

function liealgebra_so_extpowers_standard_module_struct_const(n::Int, e::Int, R::Ring) # so_n, e-th exterior power of the standard rep
    L = special_orthogonal_liealgebra(R, n)
    dimL = ngens(L)

    V = exterior_power(standard_module(L), e)
    dimV = ngens(V)
    struct_const_V = Matrix{Vector{Tuple{elem_type(R), Int}}}(undef, dimL, dimV)

    for (i, xi) in enumerate(gens(L)), (j, vj) in enumerate(gens(V))
        struct_const_V[i, j] = [(c, k) for (k, c) in enumerate(_matrix(xi * vj)) if !iszero(c)]
    end

    return struct_const_V
end


###############################################################################
#
#       generic highest weight case using GAP
#
###############################################################################

function liealgebra_gap_hightest_weight_module(dynkin::Char, n::Int, lambda::Vector{Int}, R::Ring)
    n == length(lambda) || throw(ArgumentError("length(lambda) and n have to coincide."))
    is_valid_dynkin(dynkin, n) || throw(ArgumentError("Input not allowed by GAP."))

    GAPG = GAP.Globals

    L = GAPG.SimpleLieAlgebra(GAP.julia_to_gap(string(dynkin)), n, GAPG.Rationals)
    dimL = GAPG.Dimension(L)
    basisL = GAPG.BasisVectors(GAPG.Basis(L))
    comm_table_L = GAP.gap_to_julia(GAPG.StructureConstantsTable(GAPG.Basis(L)))[1:dimL]

    struct_const_L = Matrix{Vector{Tuple{elem_type(R), Int}}}(undef, dimL, dimL)
    for i in 1:dimL, j in 1:dimL
        struct_const_L[i, j] = [(R(c), k) for (k, c) in zip(comm_table_L[i][j]...)]
    end

    V = GAPG.HighestWeightModule(L, GAP.julia_to_gap(lambda))
    dimV = GAPG.Dimension(V)
    basisV = GAPG.BasisVectors(GAPG.Basis(V))

    struct_const_V = Matrix{Vector{Tuple{elem_type(R), Int}}}(undef, dimL, dimV)
    for i in 1:dimL, j in 1:dimV
        struct_const_V[i, j] = [
            (R(c), k) for
            (k, c) in enumerate(GAP.gap_to_julia(GAPG.Coefficients(GAPG.Basis(V), basisL[i]^basisV[j]))) if !iszero(c)
        ]
    end

    symbL = ["x_$i" for i in 1:dimL]
    symbV = ["v_$i" for i in 1:dimV]

    return symbL, symbV, struct_const_L, struct_const_V
end
