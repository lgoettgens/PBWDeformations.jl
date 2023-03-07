"""
The struct representing a Lie algebra smash product.
It consists of the underlying FreeAssAlgebra with relations and some metadata.
It gets created by calling [`smash_product_lie`](@ref).
"""
@attributes mutable struct SmashProductLie{C <: RingElement}
    dimL::Int
    dimV::Int
    basisL::Vector{FreeAssAlgElem{C}}
    basisV::Vector{FreeAssAlgElem{C}}
    coeff_ring::Ring
    alg::FreeAssAlgebra{C}
    rels::QuadraticRelations{C}

    # default constructor for @attributes
    function SmashProductLie{C}(
        dimL::Int,
        dimV::Int,
        basisL::Vector{<:FreeAssAlgElem{C}},
        basisV::Vector{<:FreeAssAlgElem{C}},
        coeff_ring::Ring,
        alg::FreeAssAlgebra{C},
        rels::QuadraticRelations{C},
    ) where {C <: RingElement}
        new{C}(dimL, dimV, basisL, basisV, coeff_ring, alg, rels)
    end
end


"""
    smash_product_lie(coeff_ring::Ring, symbL::Vector{Symbol}, symbV::Vector{Symbol}, struct_const_L, struct_const_V)

Constructs the smash product over the coefficient ring `coeff_ring` using the
structure constants `struct_const_L` and `struct_const_V`, and using `symbL`
and `symbV` as symbols for the respective generators of the Lie algebra and
the module.

Returns a [`SmashProductLie`](@ref) struct and a two-part basis.
"""
function smash_product_lie(
    coeff_ring::Ring,
    symbL::Vector{Symbol},
    symbV::Vector{Symbol},
    struct_const_L::Matrix{Vector{Tuple{C, Int}}},
    struct_const_V::Matrix{Vector{Tuple{C, Int}}},
) where {C <: RingElement}
    @assert C == elem_type(coeff_ring)

    dimL = length(symbL)
    dimV = length(symbV)

    free_alg, _ = FreeAssociativeAlgebra(coeff_ring, [symbL; symbV])
    free_basisL = [gen(free_alg, i) for i in 1:dimL]
    free_basisV = [gen(free_alg, dimL + i) for i in 1:dimV]

    rels = QuadraticRelations{C}()

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

    basisL = [gen(free_alg, i) for i in 1:dimL]
    basisV = [gen(free_alg, dimL + i) for i in 1:dimV]

    return SmashProductLie{C}(dimL, dimV, basisL, basisV, coeff_ring, free_alg, rels), (basisL, basisV)
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
    struct_const_L::Matrix{Vector{Tuple{C, Int}}},
    struct_const_V::Matrix{Vector{Tuple{C, Int}}},
) where {C <: RingElem}
    return smash_product_lie(coeff_ring, map(Symbol, symbL), map(Symbol, symbV), struct_const_L, struct_const_V)
end

"""
    smash_product_lie_highest_weight(coeff_ring::Ring, dynkin::Char, n::Int, lambda::Vector{Int})

Constructs the smash product of the abstract semisimple Lie algebra given by
`dynkin` and `n` and the highest weight module with weight `lambda` over the
coefficient ring `coeff_ring`.

# Example
```jldoctest; filter = [r"(AbstractAlgebra\\.Generic\\.)?", r"(Nemo\\.)?"]
julia> smash_product_lie_highest_weight(QQ, 'A', 1, [1])
(Lie Algebra Smash Product with basis x_1, x_2, x_3, v_1, v_2 over Rational Field, (FreeAssAlgElem{fmpq}[x_1, x_2, x_3], FreeAssAlgElem{fmpq}[v_1, v_2]))
```
"""
function smash_product_lie_highest_weight(coeff_ring::Ring, dynkin::Char, n::Int, lambda::Vector{Int})
    symbL, symbV, scL, scV = liealgebra_gap_hightest_weight_module(dynkin, n, lambda, coeff_ring)

    sp, basis = smash_product_lie(coeff_ring, symbL, symbV, scL, scV)

    set_attribute!(sp, :dynkin, dynkin)
    set_attribute!(sp, :n, n)
    set_attribute!(sp, :lambda, lambda)

    return sp, basis
end

"""
    smash_product_lie_so_symmpowers_standard_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{so}_n`` and the
e-th symmetric power of the standard module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_so_symmpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # so_n, e-th symm power of standard module
    symbL = liealgebra_so_symbols(n, coeff_ring)
    scL = liealgebra_so_struct_const(n, coeff_ring)
    symbV = liealgebra_so_symmpowers_standard_module_symbols(n, e, coeff_ring)
    scV = liealgebra_so_symmpowers_standard_module_struct_const(n, e, coeff_ring)

    sp, basis = smash_product_lie(coeff_ring, symbL, symbV, scL, scV)

    set_attribute!(sp, :dynkin, n % 2 == 1 ? 'B' : 'D')
    set_attribute!(sp, :n, div(n, 2))
    set_attribute!(sp, :constructive_basis, true)
    set_attribute!(sp, :power_of_std_mod, e)

    return sp, basis
end

"""
    smash_product_lie_so_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{so}_n`` and the
e-th exterior power of the standard module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_so_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # so_n, e-th exterior power of standard module
    symbL = liealgebra_so_symbols(n, coeff_ring)
    scL = liealgebra_so_struct_const(n, coeff_ring)
    symbV = liealgebra_so_extpowers_standard_module_symbols(n, e, coeff_ring)
    scV = liealgebra_so_extpowers_standard_module_struct_const(n, e, coeff_ring)

    sp, basis = smash_product_lie(coeff_ring, symbL, symbV, scL, scV)

    set_attribute!(sp, :dynkin, n % 2 == 1 ? 'B' : 'D')
    set_attribute!(sp, :n, div(n, 2))
    set_attribute!(sp, :constructive_basis, true)
    set_attribute!(sp, :power_of_std_mod, -e)

    return sp, basis
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
    alg, _ = FreeAssociativeAlgebra(R, sp.alg.S)
    basisL = map(b -> change_base_ring(R, b, parent=alg), sp.basisL)
    basisV = map(b -> change_base_ring(R, b, parent=alg), sp.basisV)
    rels = QuadraticRelations{elem_type(R)}(k => change_base_ring(R, a, parent=alg) for (k, a) in sp.rels)

    return SmashProductLie{elem_type(R)}(sp.dimL, sp.dimV, basisL, basisV, R, alg, rels)
end
