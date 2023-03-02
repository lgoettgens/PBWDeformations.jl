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
function liealgebra_so_basis(n::Int)
    return [(b = zeros(Int, n, n); b[i, j] = 1; b[j, i] = -1; b) for i in 1:n for j in i+1:n]
end

function liealgebra_so_symbols(n::Int)
    return ["x_$(i)_$(j)" for i in 1:n for j in i+1:n]
end

function liealgebra_so_struct_const(n::Int) # so_n
    dimL = div(n * (n - 1), 2)
    basisL = liealgebra_so_basis(n)

    struct_const_L = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimL)
    for i in 1:dimL, j in 1:dimL
        struct_const_L[i, j] = [
            (c, k) for
            (k, c) in enumerate(ur_proper_triag_entries(basisL[i] * basisL[j] - basisL[j] * basisL[i])) if !iszero(c)
        ]
    end

    return struct_const_L # dimL, struct_const_L
end

function liealgebra_so_standard_module_basis(n::Int)
    return [std_basis(i, n) for i in 1:n]
end

function liealgebra_so_fundamental_module_symbols(n::Int, e::Int)
    if (n % 2 == 1 && e >= div(n, 2)) || (n % 2 == 0 && e >= div(n, 2) - 1)
        error("spin representation (not implemented) or no fundamental representation")
    end

    return liealgebra_so_extpowers_standard_module_symbols(n, e)
end

function liealgebra_so_fundamental_module_struct_const(n::Int, e::Int)
    if (n % 2 == 1 && e >= div(n, 2)) || (n % 2 == 0 && e >= div(n, 2) - 1)
        error("spin representation (not implemented) or no fundamental representation")
    end

    return liealgebra_so_extpowers_standard_module_struct_const(n, e)
end

function liealgebra_so_symmpowers_standard_module_symbols(n::Int, e::Int)
    return e == 1 ? ["v_$(i)" for i in 1:n] : ["v_$(js)" for js in Combinatorics.with_replacement_combinations(1:n, e)]
end

function liealgebra_so_symmpowers_standard_module_struct_const(n::Int, e::Int) # so_n, e-th symm power of the standard rep
    basisL = liealgebra_so_basis(n)
    basis_standardV = liealgebra_so_standard_module_basis(n)
    dimL = length(basisL)
    if e == 1
        dimV = n
        struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
        for i in 1:dimL, j in 1:dimV
            struct_const_V[i, j] = [(c, k) for (k, c) in enumerate(basisL[i] * basis_standardV[j]) if !iszero(c)]
        end
    elseif e > 1
        dimV = binomial(e + n - 1, e)
        struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
        index = Dict{Vector{Int}, Int}()
        for (j, js) in enumerate(Combinatorics.with_replacement_combinations(1:n, e))
            index[js] = j
        end
        for i in 1:dimL, js in Combinatorics.with_replacement_combinations(1:n, e)
            struct_const_V[i, index[js]] = [
                (c, index[sort([js[1:l-1]; k; js[l+1:end]])]) for l in 1:e for
                (k, c) in enumerate(basisL[i] * basis_standardV[js[l]]) if !iszero(c)
            ]
        end
    else
        error("non-negative power")
    end

    return struct_const_V # dimV, struct_const_V
end

function liealgebra_so_extpowers_standard_module_symbols(n::Int, e::Int)
    return e == 1 ? ["v_$(i)" for i in 1:n] : ["v_$(js)" for js in Combinatorics.combinations(1:n, e)]
end

function liealgebra_so_extpowers_standard_module_struct_const(n::Int, e::Int) # so_n, e-th exterior power of the standard rep
    basisL = liealgebra_so_basis(n)
    basis_standardV = liealgebra_so_standard_module_basis(n)
    dimL = length(basisL)
    if e == 1
        dimV = n
        struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
        for i in 1:dimL, j in 1:dimV
            struct_const_V[i, j] = [(c, k) for (k, c) in enumerate(basisL[i] * basis_standardV[j]) if !iszero(c)]
        end
    elseif e > 1
        dimV = binomial(n, e)
        struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
        index = Dict{Vector{Int}, Int}()
        for (j, js) in enumerate(Combinatorics.combinations(1:n, e))
            index[js] = j
        end
        for i in 1:dimL, js in Combinatorics.combinations(1:n, e)
            struct_const_V[i, index[js]] = [
                (levicivita(sortperm([js[1:l-1]; k; js[l+1:end]])) * c, index[sort([js[1:l-1]; k; js[l+1:end]])])
                for l in 1:e for (k, c) in enumerate(basisL[i] * basis_standardV[js[l]]) if
                !iszero(c) && allunique([js[1:l-1]; k; js[l+1:end]])
            ]
        end
    else
        error("non-negative power")
    end

    return struct_const_V # dimV, struct_const_V
end


###############################################################################
#
#       generic highest weight case using GAP
#
###############################################################################

function liealgebra_gap_hightest_weight_module(dynkin::Char, n::Int, lambda::Vector{Int})
    n == length(lambda) || throw(ArgumentError("length(lambda) and n have to coincide."))
    is_valid_dynkin(dynkin, n) || throw(ArgumentError("Input not allowed by GAP."))

    GAPG = GAP.Globals

    L = GAPG.SimpleLieAlgebra(GAP.julia_to_gap(string(dynkin)), n, GAPG.Rationals)
    dimL = GAPG.Dimension(L)
    basisL = GAPG.BasisVectors(GAPG.Basis(L))
    comm_table_L = GAP.gap_to_julia(GAPG.StructureConstantsTable(GAPG.Basis(L)))[1:dimL]

    struct_const_L = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimL)
    for i in 1:dimL, j in 1:dimL
        struct_const_L[i, j] = [(c, k) for (k, c) in zip(comm_table_L[i][j]...)]
    end

    V = GAPG.HighestWeightModule(L, GAP.julia_to_gap(lambda))
    dimV = GAPG.Dimension(V)
    basisV = GAPG.BasisVectors(GAPG.Basis(V))

    struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
    for i in 1:dimL, j in 1:dimV
        struct_const_V[i, j] = [
            (c, k) for
            (k, c) in enumerate(GAP.gap_to_julia(GAPG.Coefficients(GAPG.Basis(V), basisL[i]^basisV[j]))) if !iszero(c)
        ]
    end

    symbL = ["x_$i" for i in 1:dimL]
    symbV = ["v_$i" for i in 1:dimV]

    return symbL, symbV, struct_const_L, struct_const_V
end
