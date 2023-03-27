function liealgebra(R::Ring, dynkin::Tuple{Char, Int}; cached::Bool=true)
    struct_consts = liealgebra_struct_consts_gap(R, dynkin)
    L = liealgebra(R, struct_consts; cached, check=false)
    set_attribute!(L, :dynkin, dynkin)
    return L
end

function highest_weight_module(L::LieAlgebra{C}, weight::Vector{Int}; cached::Bool=true) where {C <: RingElement}
    struct_consts = liealgebra_highest_weight_module_struct_consts_gap(L, weight)
    dimV = size(struct_consts, 2)
    V = abstract_module(L, dimV, struct_consts; cached, check=false)
    set_attribute!(V, :highest_weight, weight)
    return V
end


"""
    smash_product_lie_highest_weight(coeff_ring::Ring, dynkin::Char, n::Int, lambda::Vector{Int})

Constructs the smash product of the abstract semisimple Lie algebra given by
`dynkin` and `n` and the highest weight module with weight `lambda` over the
coefficient ring `coeff_ring`.

# Example
```jldoctest; filter = [r"(AbstractAlgebra\\.Generic\\.)?", r"(Nemo\\.)?"]
julia> smash_product_lie_highest_weight(QQ, 'A', 1, [1])
Smash Product of AbstractLieAlgebra over Rational Field and AbstractModule of AbstractLieAlgebra over Rational Field
```
"""
function smash_product_lie_highest_weight(coeff_ring::Ring, dynkin::Char, n::Int, lambda::Vector{Int})
    L = liealgebra(coeff_ring, (dynkin, n))
    V = highest_weight_module(L, lambda)

    return smash_product(L, V)
end


function liealgebra_struct_consts_gap(R::Ring, dynkin::Tuple{Char, Int})
    is_valid_dynkin(dynkin...) || throw("Input not allowed by GAP.")
    R == QQ || error("Works only for QQ.")

    GAPG = GAP.Globals

    L = GAPG.SimpleLieAlgebra(GAP.julia_to_gap(string(dynkin[1])), dynkin[2], GAPG.Rationals)
    dimL = GAPG.Dimension(L)
    comm_table_L = GAP.gap_to_julia(GAPG.StructureConstantsTable(GAPG.Basis(L)))[1:dimL]

    struct_consts = Matrix{SRow{elem_type(R)}}(undef, dimL, dimL)
    for i in 1:dimL, j in 1:dimL
        struct_consts[i, j] =
            sparse_row(R, Tuple{Int, elem_type(R)}[(k, R(c)) for (k, c) in zip(comm_table_L[i][j]...)])
    end

    return struct_consts
end


function liealgebra_highest_weight_module_struct_consts_gap(
    L::LieAlgebra{C},
    weight::Vector{Int},
) where {C <: RingElement}
    R = base_ring(L)
    R == QQ || error("Works only for QQ.")
    GAPG = GAP.Globals

    gap_sc_table = [
        [
            [
                begin
                    pairs = filter(pair -> !iszero(last(pair)), collect(enumerate(_matrix(bracket(xi, xj)))))
                    [map(first, pairs), map(last, pairs)]
                end for xj in gens(L)
            ] for xi in gens(L)
        ]
        -1
        zero(base_ring(L))
    ]

    gapL = GAPG.LieAlgebraByStructureConstants(GAPG.Rationals, GAP.julia_to_gap(gap_sc_table; recursive=true))
    dimL = GAPG.Dimension(gapL)
    @assert dimL == ngens(L)
    basisL = GAPG.BasisVectors(GAPG.Basis(gapL))
    gapV = GAPG.HighestWeightModule(gapL, GAP.julia_to_gap(weight))
    dimV = GAPG.Dimension(gapV)
    basisV = GAPG.BasisVectors(GAPG.Basis(gapV))

    struct_consts = Matrix{SRow{elem_type(R)}}(undef, dimL, dimV)
    for i in 1:dimL, j in 1:dimV
        struct_consts[i, j] = sparse_row(
            R,
            Tuple{Int, elem_type(R)}[
                (k, R(c)) for
                (k, c) in enumerate(GAP.gap_to_julia(GAPG.Coefficients(GAPG.Basis(gapV), basisL[i]^basisV[j])))
            ],
        )
    end

    return struct_consts
end
