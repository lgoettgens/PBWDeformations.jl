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
