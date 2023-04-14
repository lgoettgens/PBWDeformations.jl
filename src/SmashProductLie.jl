"""
The struct representing a Lie algebra smash product.
It consists of the underlying FreeAssAlgebra with relations and some metadata.
It gets created by calling [`smash_product`](@ref).
"""
@attributes mutable struct SmashProductLie{C <: RingElement, C_lie <: RingElement}
    coeff_ring::Ring
    L::LieAlgebra{C_lie}
    V::LieAlgebraModule{C_lie}
    alg::FreeAssAlgebra{C}
    rels::QuadraticRelations{C}

    # default constructor for @attributes
    function SmashProductLie{C, C_lie}(
        coeff_ring::Ring,
        L::LieAlgebra{C_lie},
        V::LieAlgebraModule{C_lie},
        alg::FreeAssAlgebra{C},
        rels::QuadraticRelations{C},
    ) where {C <: RingElement, C_lie <: RingElement}
        new{C, C_lie}(coeff_ring, L, V, alg, rels)
    end
end

"""
    smash_product(L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElement}

Construct the smash product ``TV \\rtimes U(L)``.
"""
function smash_product(L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElement}
    L == base_lie_algebra(V) || error("Incompatible module.")
    R = base_ring(L)::parent_type(C)

    dimL = dim(L)
    dimV = dim(V)

    f_alg, _ = FreeAssociativeAlgebra(R, [symbols(L); symbols(V)])
    f_basisL = [gen(f_alg, i) for i in 1:dimL]
    f_basisV = [gen(f_alg, dimL + i) for i in 1:dimV]

    rels = QuadraticRelations{C}()

    for (i, xi) in enumerate(gens(L)), (j, xj) in enumerate(gens(L))
        commutator =
            sum(c * f_basisL[k] for (k, c) in enumerate(_matrix(bracket(xi, xj))) if !iszero(c); init=zero(f_alg))
        rels[(i, j)] = f_basisL[j] * f_basisL[i] + commutator

    end

    for (i, xi) in enumerate(gens(L)), (j, vj) in enumerate(gens(V))
        commutator = sum(c * f_basisV[k] for (k, c) in enumerate(_matrix(xi * vj)) if !iszero(c); init=zero(f_alg))
        rels[(i, dimL + j)] = f_basisV[j] * f_basisL[i] + commutator
        rels[(dimL + j, i)] = f_basisL[i] * f_basisV[j] - commutator
    end

    sp = SmashProductLie{C, C}(R, L, V, f_alg, rels)

    return sp
end


ngens(sp::SmashProductLie{C}) where {C <: RingElement} = length(gens(sp.alg)) # ngens(sp.alg), see https://github.com/Nemocas/AbstractAlgebra.jl/pull/1295
function ngens(sp::SmashProductLie{C}, part::Symbol) where {C <: RingElement}
    part == :L && return dim(sp.L)
    part == :V && return dim(sp.V)
    error("Invalid part.")
end

gens(sp::SmashProductLie{C}) where {C <: RingElement} = gens(sp.alg)
function gens(sp::SmashProductLie{C}, part::Symbol) where {C}
    part == :L && return [gen(sp.alg, i) for i in 1:dim(sp.L)]
    part == :V && return [gen(sp.alg, i + dim(sp.L)) for i in 1:dim(sp.V)]
    error("Invalid part.")
end

gen(sp::SmashProductLie{C}, i::Int) where {C <: RingElement} = gen(sp.alg, i)
function gen(sp::SmashProductLie{C}, i::Int, part::Symbol) where {C <: RingElement}
    1 <= i <= ngens(sp, part) || error("Invalid generator index.")
    part == :L && return gen(sp.alg, i)
    part == :V && return gen(sp.alg, i + dim(sp.L))
    error("Invalid part.")
end

function show(io::IO, sp::SmashProductLie{C, C_lie}) where {C <: RingElement, C_lie <: RingElement}
    print(io, "Smash Product")
    if C_lie != C
        print(io, " over ")
        print(IOContext(io, :compact => true), base_ring(sp.alg))
    end
    print(io, " of ")
    print(IOContext(io, :compact => true), sp.L)
    print(io, " and ")
    print(IOContext(io, :compact => true), sp.V)
end


function change_base_ring(R::Ring, sp::SmashProductLie{C, C_lie}) where {C <: RingElement, C_lie <: RingElement}
    alg, _ = FreeAssociativeAlgebra(R, symbols(sp.alg))
    rels = QuadraticRelations{elem_type(R)}(k => change_base_ring(R, a, parent=alg) for (k, a) in sp.rels)

    return SmashProductLie{elem_type(R), C_lie}(R, sp.L, sp.V, alg, rels)
end
