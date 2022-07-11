"""
The struct representing a lie algebra smash product.
It consists of the underlying QuadraticQuoAlgebra and some metadata.
"""
mutable struct SmashProductLie{C <: RingElement}
    dimL::Int64
    dimV::Int64
    basisL::Vector{QuadraticQuoAlgebraElem{C}}
    basisV::Vector{QuadraticQuoAlgebraElem{C}}
    coeff_ring::Ring
    alg::QuadraticQuoAlgebra{C}
    # dynkin :: Char
    # n :: Int64
    # lambda :: Vector{Int64}
    # matrixRepL :: Vector{Matrix{Int64}}
end


"""
    smash_product_lie(coeff_ring::Ring, symbL::Vector{Symbol}, symbV::Vector{Symbol}, struct_const_L, struct_const_V)

Constructs the smash product over the coefficient ring `coeff_ring` using the
structure constants `struct_const_L` and `struct_const_V`, and using `symbL`
and `symbV` as symbols for the respective generators of the lie algebra and
the module.
"""
function smash_product_lie(
    coeff_ring::Ring,
    symbL::Vector{Symbol},
    symbV::Vector{Symbol},
    struct_const_L::Matrix{Vector{Tuple{Int, Int}}},
    struct_const_V::Matrix{Vector{Tuple{Int, Int}}},
)
    C = elem_type(coeff_ring)

    dimL = length(symbL)
    dimV = length(symbV)

    free_alg, _ = free_algebra(coeff_ring, [symbL; symbV])
    free_basisL = [gen(free_alg, i) for i in 1:dimL]
    free_basisV = [gen(free_alg, dimL + i) for i in 1:dimV]

    rels = Dict{Tuple{Int, Int}, FreeAlgebraElem{C}}()

    for i in 1:dimL, j in 1:dimL
        rels[(i, j)] =
            free_basisL[j] * free_basisL[i] +
            sum(c * free_basisL[k] for (c, k) in struct_const_L[i, j]; init=zero(free_alg))
    end

    for i in 1:dimL, j in 1:dimV
        rels[(i, dimL + j)] =
            free_basisV[j] * free_basisL[i] +
            sum(c * free_basisV[k] for (c, k) in struct_const_V[i, j]; init=zero(free_alg))
        rels[(dimL + j, i)] =
            free_basisL[i] * free_basisV[j] -
            sum(c * free_basisV[k] for (c, k) in struct_const_V[i, j]; init=zero(free_alg))
    end

    alg, _ = quadratic_quo_algebra(free_alg, rels)
    basisL = [gen(alg, i) for i in 1:dimL]
    basisV = [gen(alg, dimL + i) for i in 1:dimV]

    return SmashProductLie{C}(dimL, dimV, basisL, basisV, coeff_ring, alg), (basisL, basisV)
end

"""
    smash_product_lie(coeff_ring::Ring, symbL::Vector{String}, symbV::Vector{String}, struct_const_L, struct_const_V)

The same as the other method with structure constants, but takes strings
instead of symbols to name the generators.
"""
function smash_product_lie(
    coeff_ring::Ring,
    symbL::Vector{String},
    symbV::Vector{String},
    struct_const_L::Matrix{Vector{Tuple{Int, Int}}},
    struct_const_V::Matrix{Vector{Tuple{Int, Int}}},
)
    return smash_product_lie(coeff_ring, map(Symbol, symbL), map(Symbol, symbV), struct_const_L, struct_const_V)
end

"""
    smash_product_lie(coeff_ring::Ring, dynkin::Char, n::Int, lambda::Vector{Int})

Constructs the smash product of the abstract semisimple lie algebra given by
`dynkin` and `n` and the highest weight module with weight `lambda` over the
coefficient ring `coeff_ring`.

# Example
```jldoctest
julia> smash_product_lie(QQ, 'A', 1, [1])
(Lie Algebra Smash Product with basis x_1, x_2, x_3, v_1, v_2 over Rational Field, (QuadraticQuoAlgebraElem{fmpq}[x_1, x_2, x_3], QuadraticQuoAlgebraElem{fmpq}[v_1, v_2]))
```
"""
function smash_product_lie(coeff_ring::Ring, dynkin::Char, n::Int, lambda::Vector{Int})
    dimL, dimV, struct_const_L, struct_const_V = smash_product_struct_const_from_gap(dynkin, n, lambda)

    symbL = ["x_$i" for i in 1:dimL]
    symbV = ["v_$i" for i in 1:dimV]

    return smash_product_lie(coeff_ring, symbL, symbV, struct_const_L, struct_const_V)
end

function smash_product_struct_const_from_gap(dynkin::Char, n::Int, lambda::Vector{Int})
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

    return dimL, dimV, struct_const_L, struct_const_V
end

"""
    smash_product_lie_so(coeff_ring :: Ring, n :: Int, lambda :: Vector{Int})

Constructs the smash product of the semisimple lie algebra `so_n` and the highest weight module with weight `lambda` over the
coefficient ring `coeff_ring`. For `n` odd, `lambda` is a B-type highest weight; for `n` even, `lambda` is a D-type highest weight.

# Example
```jldoctest
julia> smash_product_lie_so(QQ, 5, [1,0])
(Lie Algebra Smash Product with basis x_1_2, x_1_3, x_1_4, ..., x_4_5, v_1, v_2, v_3, ..., v_5 over Rational Field, (QuadraticQuoAlgebraElem{fmpq}[x_1_2, x_1_3, x_1_4, x_1_5, x_2_3, x_2_4, x_2_5, x_3_4, x_3_5, x_4_5], QuadraticQuoAlgebraElem{fmpq}[v_1, v_2, v_3, v_4, v_5]))
```

!!! warning
    Only fundamental non-spin representations are currently supported.
"""
function smash_product_lie_so(coeff_ring::Ring, n::Int, lambda::Vector{Int}) # for odd n this is a B type highest weight, for even n D type
    _, _, struct_const_L, struct_const_V, symbL, symbV = smash_product_struct_const_so(n, lambda)

    return smash_product_lie(coeff_ring, symbL, symbV, struct_const_L, struct_const_V)
end

function smash_product_struct_const_so(n::Int, lambda::Vector{Int}) # for odd n this is a B type highest weight, for even n D type
    @assert div(n, 2) == length(lambda)
    lenghtlambda = div(n, 2)

    ur_triag(M) = vcat([M[i, i+1:end] for i in 1:size(M, 1)]...)
    std_basis(i, n) = [i == j ? 1 : 0 for j in 1:n]

    dimL = div(n * (n - 1), 2)
    basisL = [(b = zeros(Int, n, n); b[i, j] = 1; b[j, i] = -1; b) for i in 1:n for j in i+1:n]
    symbL = ["x_$(i)_$(j)" for i in 1:n for j in i+1:n]

    struct_const_L = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimL)
    for i in 1:dimL, j in 1:dimL
        struct_const_L[i, j] =
            [(c, k) for (k, c) in enumerate(ur_triag(basisL[i] * basisL[j] - basisL[j] * basisL[i])) if !iszero(c)]
    end

    if in(lambda, [std_basis(i, lenghtlambda) for i in 1:lenghtlambda])
        lambda_std_i = findfirst(==(lambda), [std_basis(i, lenghtlambda) for i in 1:lenghtlambda])

        if (n % 2 == 1 && lambda_std_i <= lenghtlambda - 1) || (n % 2 == 0 && lambda_std_i <= lenghtlambda - 2) # exterior product of defining representation
            if lambda_std_i == 1
                dimV = n
                symbV = ["v_$(i)" for i in 1:dimV]
                struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
                for i in 1:dimL, j in 1:dimV
                    struct_const_V[i, j] = [(c, k) for (k, c) in enumerate(basisL[i] * std_basis(j, n)) if !iszero(c)]
                end
            else
                dimV = binomial(n, lambda_std_i)
                symbV = ["v_$(js)" for js in Combinatorics.combinations(1:n, lambda_std_i)] # TODO
                struct_const_V = Matrix{Vector{Tuple{Int, Int}}}(undef, dimL, dimV)
                index = Dict{Vector{Int}, Int}()
                for (j, js) in enumerate(Combinatorics.combinations(1:n, lambda_std_i))
                    index[js] = j
                end
                for i in 1:dimL, js in Combinatorics.combinations(1:n, lambda_std_i)
                    struct_const_V[i, index[js]] = [
                        (
                            levicivita(sortperm([js[1:l-1]; k; js[l+1:end]])) * c,
                            index[sort([js[1:l-1]; k; js[l+1:end]])],
                        ) for l in 1:lambda_std_i for (k, c) in enumerate(basisL[i] * std_basis(js[l], n)) if
                        !iszero(c) && allunique([js[1:l-1]; k; js[l+1:end]])
                    ]
                end
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

function gens(sp::SmashProductLie{C}) where {C <: RingElement}
    return [gen(sp.alg, i) for i in 1:sp.dimL], [gen(sp.alg, i + sp.dimL) for i in 1:sp.dimV]
end


function show(io::IO, sp::SmashProductLie)
    local max_gens = 4 # largest number of generators to print
    print(io, "Lie Algebra Smash Product with basis ")
    for i in 1:min(sp.dimL - 1, max_gens - 1)
        print(io, string(sp.alg.S[i]), ", ")
    end
    if sp.dimL > max_gens
        print(io, "..., ")
    end
    print(io, string(sp.alg.S[sp.dimL]) * ", ")
    for i in 1:min(sp.dimV - 1, max_gens - 1)
        print(io, string(sp.alg.S[sp.dimL+i]), ", ")
    end
    if sp.dimV > max_gens
        print(io, "..., ")
    end
    print(io, string(sp.alg.S[sp.dimL+sp.dimV]))
    print(io, " over ")
    print(IOContext(io, :compact => true), sp.coeff_ring)
end


function change_base_ring(R::Ring, sp::SmashProductLie{C}) where {C <: RingElement}
    alg = change_base_ring(R, sp.alg)
    basisL = [gen(alg, i) for i in 1:sp.dimL]
    basisV = [gen(alg, sp.dimL + i) for i in 1:sp.dimV]

    return SmashProductLie{elem_type(R)}(sp.dimL, sp.dimV, basisL, basisV, R, alg)
end
