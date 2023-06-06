"""
The struct representing a Lie algebra smash product.
It consists of the underlying FreeAssAlgebra with relations and some metadata.
It gets created by calling [`smash_product`](@ref).
"""
@attributes mutable struct SmashProductLie{C <: RingElem, CL <: RingElem} <: NCRing
    coeff_ring::Ring
    L::LieAlgebra{CL}
    V::LieAlgebraModule{CL}
    alg::FreeAssAlgebra{C}
    rels::Matrix{Union{Nothing, FreeAssAlgElem{C}}}

    # default constructor for @attributes
    function SmashProductLie{C, CL}(
        coeff_ring::Ring,
        L::LieAlgebra{CL},
        V::LieAlgebraModule{CL},
        alg::FreeAssAlgebra{C},
        rels::Matrix{Union{Nothing, FreeAssAlgElem{C}}},
    ) where {C <: RingElem, CL <: RingElem}
        new{C, CL}(coeff_ring, L, V, alg, rels)
    end
end

mutable struct SmashProductLieElem{C <: RingElem, CL <: RingElem} <: NCRingElem
    p::SmashProductLie{C, CL}   # parent
    alg_elem::FreeAssAlgElem{C}
    simplified::Bool

    function SmashProductLieElem(
        p::SmashProductLie{C, CL},
        alg_elem::FreeAssAlgElem{C};
        simplified::Bool=false,
    ) where {C <: RingElem, CL <: RingElem}
        @req underlying_algebra(p) === parent(alg_elem) "Incompatible algebras."
        return new{C, CL}(p, alg_elem, simplified)
    end
end

###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{SmashProductLieElem{C, CL}}) where {C <: RingElem, CL <: RingElem} = SmashProductLie{C, CL}

elem_type(::Type{SmashProductLie{C, CL}}) where {C <: RingElem, CL <: RingElem} = SmashProductLieElem{C, CL}

parent(e::SmashProductLieElem) = e.p

base_ring(Sp::SmashProductLie) = Sp.coeff_ring

base_ring(e::SmashProductLieElem) = base_ring(parent(e))

coefficient_ring(Sp::SmashProductLie) = Sp.coeff_ring

lie_algebra(Sp::SmashProductLie) = Sp.L

lie_module(Sp::SmashProductLie) = Sp.V

underlying_algebra(Sp::SmashProductLie) = Sp.alg

ngens(Sp::SmashProductLie) = ngens(underlying_algebra(Sp))
function ngens(Sp::SmashProductLie, part::Symbol)
    part == :L && return dim(lie_algebra(Sp))
    part == :V && return dim(lie_module(Sp))
    error("Invalid part.")
end

gens(Sp::SmashProductLie) = map(Sp, gens(underlying_algebra(Sp)))
function gens(Sp::SmashProductLie, part::Symbol)
    part == :L && return [Sp(gen(underlying_algebra(Sp), i)) for i in 1:dim(lie_algebra(Sp))]
    part == :V && return [Sp(gen(underlying_algebra(Sp), i + dim(lie_algebra(Sp)))) for i in 1:dim(lie_module(Sp))]
    error("Invalid part.")
end

gen(Sp::SmashProductLie, i::Int) = Sp(gen(underlying_algebra(Sp), i))
function gen(Sp::SmashProductLie, i::Int, part::Symbol)
    @req 1 <= i <= ngens(Sp, part) "Invalid generator index."
    part == :L && return Sp(gen(underlying_algebra(Sp), i))
    part == :V && return Sp(gen(underlying_algebra(Sp), i + dim(lie_algebra(Sp))))
    error("Invalid part.")
end

function zero(Sp::SmashProductLie)
    return Sp(zero(underlying_algebra(Sp)))
end

function iszero(e::SmashProductLieElem)
    return iszero(simplify(e).alg_elem)
end

function one(Sp::SmashProductLie)
    return Sp(one(underlying_algebra(Sp)))
end

function isone(e::SmashProductLieElem)
    return isone(simplify(e).alg_elem)
end

function Base.deepcopy_internal(e::SmashProductLieElem, dict::IdDict)
    return SmashProductLieElem(parent(e), deepcopy_internal(e.alg_elem, dict); simplified=e.simplified)
end

function check_parent(e1::SmashProductLieElem{C}, e2::SmashProductLieElem{C}) where {C <: RingElem}
    parent(e1) != parent(e2) && error("Incompatible smash products.")
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, Sp::SmashProductLie{C, CL}) where {C <: RingElem, CL <: RingElem}
    print(io, "Smash Product")
    if CL != C
        print(io, " over ")
        print(IOContext(io, :supercompact => true), base_ring(underlying_algebra(Sp)))
    end
    print(io, " of ")
    print(IOContext(io, :compact => true), lie_algebra(Sp))
    print(io, " and ")
    print(IOContext(io, :compact => true), lie_module(Sp))
end


function show(io::IO, e::SmashProductLieElem)
    show(io, e.alg_elem)
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (Sp::SmashProductLie)()
    return zero(Sp)
end

function (Sp::SmashProductLie)(e::Union{RingElement, NCRingElem})
    return Sp(underlying_algebra(Sp)(e))
end

function (Sp::SmashProductLie{C, CL})(e::FreeAssAlgElem{C}) where {C <: RingElem, CL <: RingElem}
    if underlying_algebra(Sp) !== parent(e)
        e = underlying_algebra(Sp)(e)
    end
    return SmashProductLieElem(Sp, e)
end

function (Sp::SmashProductLie{C, CL})(e::SmashProductLieElem{C, CL}) where {C <: RingElem, CL <: RingElem}
    @req parent(e) == Sp "Incompatible smash products."
    return e
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(e::SmashProductLieElem)
    return parent(e)(-e.alg_elem)
end

function Base.:+(e1::SmashProductLieElem{C}, e2::SmashProductLieElem{C}) where {C <: RingElem}
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem + e2.alg_elem)
end

function Base.:-(e1::SmashProductLieElem{C}, e2::SmashProductLieElem{C}) where {C <: RingElem}
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem - e2.alg_elem)
end

function Base.:*(e1::SmashProductLieElem{C}, e2::SmashProductLieElem{C}) where {C <: RingElem}
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem * e2.alg_elem)
end

function Base.:*(e::SmashProductLieElem{C}, c::C) where {C <: RingElem}
    base_ring(e) != parent(c) && error("Incompatible rings.")
    return parent(e)(e.alg_elem * c)
end

function Base.:*(e::SmashProductLieElem{C}, c::U) where {C <: RingElem, U <: Union{Rational, Integer}}
    return parent(e)(e.alg_elem * c)
end

function Base.:*(c::C, e::SmashProductLieElem{C}) where {C <: RingElem}
    base_ring(e) != parent(c) && error("Incompatible rings.")
    return parent(e)(c * e.alg_elem)
end

function Base.:*(c::U, e::SmashProductLieElem{C}) where {C <: RingElem, U <: Union{Rational, Integer}}
    return parent(e)(c * e.alg_elem)
end

function Base.:^(e::SmashProductLieElem, n::Int)
    return parent(e)(e.alg_elem^n)
end

function comm(e1::SmashProductLieElem{C}, e2::SmashProductLieElem{C}) where {C <: RingElem}
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem * e2.alg_elem - e2.alg_elem * e1.alg_elem)
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(e1::SmashProductLieElem, e2::SmashProductLieElem)
    return parent(e1) === parent(e2) && simplify(e1).alg_elem == simplify(e2).alg_elem
end

function Base.hash(e::SmashProductLieElem, h::UInt)
    e = simplify(e)
    b = 0xdcc11ff793ca4ada % UInt
    h = hash(parent(e), h)
    h = hash(e.alg_elem, h)
    return xor(h, b)
end

###############################################################################
#
#   Simplification
#
###############################################################################

function simplify(e::SmashProductLieElem)
    e.simplified && return e
    e.alg_elem = _normal_form(e.alg_elem, parent(e).rels)
    e.simplified = true
    return e
end

function _normal_form(a::FreeAssAlgElem{C}, rels::Matrix{Union{Nothing, FreeAssAlgElem{C}}}) where {C <: RingElem}
    todo = deepcopy(a)
    result = zero(parent(todo))
    CR = base_ring(a)
    A = parent(todo)
    while todo.length > 0
        c = leading_coefficient(todo)
        exp = leading_exponent_word(todo)
        t = leading_term(todo)
        todo -= t

        changed = false
        for i in 1:length(exp)-1
            if exp[i] > exp[i+1] && !isnothing(rels[exp[i], exp[i+1]])
                changed = true
                todo += A([c], [exp[1:i-1]]) * rels[exp[i], exp[i+1]] * A([one(CR)], [exp[i+2:end]])
                break
            end
        end
        if !changed
            result += t
        end
    end
    return result
end

###############################################################################
#
#   Constructor
#
###############################################################################

"""
    smash_product(L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElem}

Construct the smash product ``TV \\rtimes U(L)``.
"""
function smash_product(L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElem}
    return smash_product(base_ring(L), L, V)
end

"""
    smash_product(R::Ring, L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElem}

Construct the smash product ``TV \\rtimes U(L)`` and extend the coefficients to `R`.
"""
function smash_product(R::Ring, L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElem}
    @req L == base_lie_algebra(V) "Incompatible module."
    try
        R(zero(base_ring(L)))
    catch e
        if e isa MethodError
            error("Coefficient extension from $(base_ring(L)) to $(R) not possible.")
        else
            rethrow(e)
        end
    end

    dimL = dim(L)
    dimV = dim(V)

    f_alg, _ = free_associative_algebra(
        R,
        [symbols(L); is_standard_module(V) ? symbols(V) : (x -> Symbol("($x)")).(symbols(V))],
    )
    f_basisL = [gen(f_alg, i) for i in 1:dimL]
    f_basisV = [gen(f_alg, dimL + i) for i in 1:dimV]

    rels = Matrix{Union{Nothing, FreeAssAlgElem{elem_type(R)}}}(nothing, dimL + dimV, dimL + dimV)

    for (i, xi) in enumerate(basis(L)), (j, xj) in enumerate(basis(L))
        commutator = sum(R(c) * f_basisL[k] for (k, c) in enumerate(_matrix(xi * xj)) if !iszero(c); init=zero(f_alg))
        rels[i, j] = f_basisL[j] * f_basisL[i] + commutator

    end

    for (i, xi) in enumerate(basis(L)), (j, vj) in enumerate(basis(V))
        commutator = sum(R(c) * f_basisV[k] for (k, c) in enumerate(_matrix(xi * vj)) if !iszero(c); init=zero(f_alg))
        rels[i, dimL+j] = f_basisV[j] * f_basisL[i] + commutator
        rels[dimL+j, i] = f_basisL[i] * f_basisV[j] - commutator
    end

    Sp = SmashProductLie{elem_type(R), C}(R, L, V, f_alg, rels)

    return Sp
end
