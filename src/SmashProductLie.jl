"""
A struct containing additional information about a Lie algebra smash product.
Every field can be set to nothing if it is unknown.

| Name | Type | Description |
|:---- |:---- |:------------|
| `dynkin` | `Char?` | the family of the dynkin type of the Lie algebra |
| `n` | `Int?`| the `n` of the dynkin type of the Lie algebra |
| `lambda` | `Vector{Int}?`| highest weight vector of the module, only existing if the module is simple |
| `constructive_basis` | `Bool`| true if the used basis for the structure constants is known in terms of matrices |
| `power_of_std_mod` | `Int?`| if the module is a power of the standard module, positive = symmetric power, negative = exterior power |
"""
struct SmashProductLieInfo
    dynkin::Union{Nothing, Char}
    n::Union{Nothing, Int}
    lambda::Union{Nothing, Vector{Int}}     # highest weight vector
    constructive_basis::Bool
    power_of_std_mod::Union{Nothing, Int}   # positive = symmetric power, negative = exterior power
    function SmashProductLieInfo(;
        dynkin=nothing,
        n=nothing,
        lambda=nothing,
        constructive_basis=false,
        power_of_std_mod=nothing,
    )
        new(dynkin, n, lambda, constructive_basis, power_of_std_mod)
    end
end


"""
The struct representing a Lie algebra smash product.
It consists of the underlying FreeAlgebra with relations and some metadata.
It gets created by calling [`smash_product_lie`](@ref).
"""
mutable struct SmashProductLie{C <: RingElement}
    dimL::Int
    dimV::Int
    basisL::Vector{FreeAlgebraElem{C}}
    basisV::Vector{FreeAlgebraElem{C}}
    coeff_ring::Ring
    alg::FreeAlgebra{C}
    rels::QuadraticRelations{C}
    info::SmashProductLieInfo
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
    struct_const_L::Matrix{Vector{Tuple{Int, Int}}},
    struct_const_V::Matrix{Vector{Tuple{Int, Int}}},
    info=SmashProductLieInfo()::SmashProductLieInfo,
)
    C = elem_type(coeff_ring)

    dimL = length(symbL)
    dimV = length(symbV)

    free_alg, _ = free_algebra(coeff_ring, [symbL; symbV])
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

    return SmashProductLie{C}(dimL, dimV, basisL, basisV, coeff_ring, free_alg, rels, info), (basisL, basisV)
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
    info=SmashProductLieInfo()::SmashProductLieInfo,
)
    return smash_product_lie(coeff_ring, map(Symbol, symbL), map(Symbol, symbV), struct_const_L, struct_const_V, info)
end

"""
    smash_product_lie_highest_weight(coeff_ring::Ring, dynkin::Char, n::Int, lambda::Vector{Int})

Constructs the smash product of the abstract semisimple Lie algebra given by
`dynkin` and `n` and the highest weight module with weight `lambda` over the
coefficient ring `coeff_ring`.

# Example
```jldoctest
julia> smash_product_lie_highest_weight(QQ, 'A', 1, [1])
(Lie Algebra Smash Product with basis x_1, x_2, x_3, v_1, v_2 over Rational Field, (FreeAlgebraElem{fmpq}[x_1, x_2, x_3], FreeAlgebraElem{fmpq}[v_1, v_2]))
```
"""
function smash_product_lie_highest_weight(coeff_ring::Ring, dynkin::Char, n::Int, lambda::Vector{Int})
    symbL, symbV, scL, scV = liealgebra_gap_hightest_weight_module(dynkin, n, lambda)

    info = SmashProductLieInfo(dynkin=dynkin, n=n, lambda=lambda)

    return smash_product_lie(coeff_ring, symbL, symbV, scL, scV, info)
end

"""
    smash_product_lie_so_fundamental_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{so}_n`` and the
e-th fundamental module over the coefficient ring `coeff_ring`.
"""
function smash_product_lie_so_fundamental_module(coeff_ring::Ring, n::Int, e::Int) # so_n, e-th fundamental module (spin reps not implemented)
    symbL = liealgebra_so_symbols(n)
    scL = liealgebra_so_struct_const(n)
    symbV = liealgebra_so_fundamental_module_symbols(n, e)
    scV = liealgebra_so_fundamental_module_struct_const(n, e)

    info = SmashProductLieInfo(dynkin=(n % 2 == 1 ? 'B' : 'D'), n=div(n, 2), constructive_basis=true)

    return smash_product_lie(coeff_ring, symbL, symbV, scL, scV, info)
end

"""
    smash_product_lie_so_symmpowers_fundamental_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{so}_n`` and the
e-th symmetric power of the fundamental module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_so_symmpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # so_n, e-th symm power of standard module
    symbL = liealgebra_so_symbols(n)
    scL = liealgebra_so_struct_const(n)
    symbV = liealgebra_so_symmpowers_standard_module_symbols(n, e)
    scV = liealgebra_so_symmpowers_standard_module_struct_const(n, e)

    info =
        SmashProductLieInfo(dynkin=(n % 2 == 1 ? 'B' : 'D'), n=div(n, 2), constructive_basis=true, power_of_std_mod=e)

    return smash_product_lie(coeff_ring, symbL, symbV, scL, scV, info)
end

"""
    smash_product_lie_so_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{so}_n`` and the
e-th exterior power of the fundamental module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_so_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # so_n, e-th exterior power of standard module
    symbL = liealgebra_so_symbols(n)
    scL = liealgebra_so_struct_const(n)
    symbV = liealgebra_so_extpowers_standard_module_symbols(n, e)
    scV = liealgebra_so_extpowers_standard_module_struct_const(n, e)

    info =
        SmashProductLieInfo(dynkin=(n % 2 == 1 ? 'B' : 'D'), n=div(n, 2), constructive_basis=true, power_of_std_mod=-e)

    return smash_product_lie(coeff_ring, symbL, symbV, scL, scV, info)
end

"""
    smash_product_lie_sp_symmpowers_standard_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{sp}_{2n}`` and the
e-th symmetric power of the standard module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_sp_symmpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # sp_2n, e-th symm power of standard module
    symbL = liealgebra_sp_symbols(n)
    scL = liealgebra_sp_struct_const(n)
    symbV = liealgebra_sp_symmpowers_standard_module_symbols(n, e)
    scV = liealgebra_sp_symmpowers_standard_module_struct_const(n, e)

    info = SmashProductLieInfo(dynkin='C', n=n, constructive_basis=true, power_of_std_mod=e)

    return smash_product_lie(coeff_ring, symbL, symbV, scL, scV, info)
end

"""
    smash_product_lie_sp_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int)

Constructs the smash product of the Lie algebra ``\\mathfrak{sp}_{2n}`` and the
e-th exterior power of the fundamental module over the
coefficient ring `coeff_ring`.
"""
function smash_product_lie_sp_extpowers_standard_module(coeff_ring::Ring, n::Int, e::Int) # sp_2n, e-th exterior power of standard module
    symbL = liealgebra_sp_symbols(n)
    scL = liealgebra_sp_struct_const(n)
    symbV = liealgebra_sp_extpowers_standard_module_symbols(n, e)
    scV = liealgebra_sp_extpowers_standard_module_struct_const(n, e)

    info = SmashProductLieInfo(dynkin='C', n=n, constructive_basis=true, power_of_std_mod=-e)

    return smash_product_lie(coeff_ring, symbL, symbV, scL, scV, info)
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
    rels = QuadraticRelations{elem_type(R)}(k => alg(a) for (k, a) in sp.rels)

    return SmashProductLie{elem_type(R)}(sp.dimL, sp.dimV, basisL, basisV, R, alg, rels, sp.info)
end
