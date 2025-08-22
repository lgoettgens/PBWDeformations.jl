###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(
    ::Type{SmashProductLieElem{C, LieC, LieT}},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = SmashProductLie{C, LieC, LieT}

elem_type(
    ::Type{SmashProductLie{C, LieC, LieT}},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = SmashProductLieElem{C, LieC, LieT}

parent(e::SmashProductLieElem) = e.p

coefficient_ring(
    Sp::SmashProductLie{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = Sp.coeff_ring::parent_type(C)

coefficient_ring(e::SmashProductLieElem) = coefficient_ring(parent(e))

base_lie_algebra(
    Sp::SmashProductLie{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = Sp.L::parent_type(LieT)

base_module(Sp::SmashProductLie{C, LieC, LieT}) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} =
    Sp.V

underlying_algebra(
    Sp::SmashProductLie{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = Sp.alg::free_associative_algebra_type(C)

data(e::SmashProductLieElem{C, LieC, LieT}) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} =
    e.alg_elem::elem_type(free_associative_algebra_type(C))

ngens(Sp::SmashProductLie) = ngens(underlying_algebra(Sp))
function ngens(Sp::SmashProductLie, part::Symbol)
    part == :L && return dim(base_lie_algebra(Sp))
    part == :V && return dim(base_module(Sp))
    error("Invalid part.")
end

gens(Sp::SmashProductLie) = map(Sp, gens(underlying_algebra(Sp)))
function gens(Sp::SmashProductLie, part::Symbol)
    part == :L && return [Sp(gen(underlying_algebra(Sp), i)) for i in 1:dim(base_lie_algebra(Sp))]
    part == :V &&
        return [Sp(gen(underlying_algebra(Sp), i + dim(base_lie_algebra(Sp)))) for i in 1:dim(base_module(Sp))]
    error("Invalid part.")
end

gen(Sp::SmashProductLie, i::Int) = Sp(gen(underlying_algebra(Sp), i))
function gen(Sp::SmashProductLie, i::Int, part::Symbol)
    @req 1 <= i <= ngens(Sp, part) "Invalid generator index."
    part == :L && return Sp(gen(underlying_algebra(Sp), i))
    part == :V && return Sp(gen(underlying_algebra(Sp), i + dim(base_lie_algebra(Sp))))
    error("Invalid part.")
end

function zero(Sp::SmashProductLie)
    return Sp(zero(underlying_algebra(Sp)))
end

function iszero(e::SmashProductLieElem)
    return iszero(data(simplify(e)))
end

function one(Sp::SmashProductLie)
    return Sp(one(underlying_algebra(Sp)))
end

function isone(e::SmashProductLieElem)
    return isone(data(simplify(e)))
end

function Base.deepcopy_internal(e::SmashProductLieElem, dict::IdDict)
    return SmashProductLieElem(parent(e), deepcopy_internal(data(e), dict); simplified=e.simplified)
end

function check_parent(
    e1::SmashProductLieElem{C, LieC, LieT},
    e2::SmashProductLieElem{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    parent(e1) != parent(e2) && error("Incompatible smash products.")
end

function AbstractAlgebra.promote_rule(
    T1::Type{<:SmashProductLieElem{C1}},
    ::Type{C2},
) where {C1 <: RingElem, C2 <: RingElem}
    AbstractAlgebra.promote_rule(C1, C2) == C1 ? T1 : Union{}
end

function change_base_ring(R::Ring, e::SmashProductLieElem{C}; parent::SmashProductLie=smash_product(R, base_lie_algebra(parent(e)), base_module(parent(e)))) where C
    return parent(change_base_ring(R, data(e); parent=underlying_algebra(parent)))
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, Sp::SmashProductLie{C, LieC}) where {C <: RingElem, LieC <: FieldElem}
    print(io, "Smash Product")
    if LieC != C
        print(io, " over ")
        print(IOContext(io, :supercompact => true), coefficient_ring(underlying_algebra(Sp)))
    end
    print(io, " of ")
    print(IOContext(io, :compact => true), base_lie_algebra(Sp))
    print(io, " and ")
    print(IOContext(io, :compact => true), base_module(Sp))
end


function show(io::IO, e::SmashProductLieElem)
    show(io, data(e))
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

function (Sp::SmashProductLie{C})(e::FreeAssociativeAlgebraElem{C}) where {C <: RingElem}
    if underlying_algebra(Sp) !== parent(e)
        e = underlying_algebra(Sp)(e)
    end
    return SmashProductLieElem(Sp, e)
end

function (Sp::SmashProductLie{C, LieC, LieT})(
    e::SmashProductLieElem{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req parent(e) == Sp "Incompatible smash products."
    return e
end

function (Sp::SmashProductLie{C, LieC, LieT})(
    x::LieT,
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req parent(x) == base_lie_algebra(Sp) "Incompatible smash products."
    return sum(c * b for (c, b) in zip(coefficients(x), gens(Sp, :L)))
end

function (Sp::SmashProductLie{C, LieC, LieT})(
    v::LieAlgebraModuleElem{LieC},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req parent(v) == base_module(Sp) "Incompatible smash products."
    return sum(c * b for (c, b) in zip(coefficients(v), gens(Sp, :V)))
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(e::SmashProductLieElem)
    return parent(e)(-data(e))
end

function Base.:+(e1::SmashProductLieElem, e2::SmashProductLieElem)
    check_parent(e1, e2)
    return parent(e1)(data(e1) + data(e2))
end

function Base.:-(e1::SmashProductLieElem, e2::SmashProductLieElem)
    check_parent(e1, e2)
    return parent(e1)(data(e1) - data(e2))
end

function Base.:*(e1::SmashProductLieElem, e2::SmashProductLieElem)
    check_parent(e1, e2)
    return parent(e1)(data(e1) * data(e2))
end

function Base.:*(e::SmashProductLieElem{C}, c::C) where {C <: RingElem}
    coefficient_ring(e) != parent(c) && error("Incompatible rings.")
    return parent(e)(data(e) * c)
end

function Base.:*(e::SmashProductLieElem, c::U) where {U <: Union{Rational, IntegerUnion}}
    return parent(e)(data(e) * c)
end

function Base.:*(e::SmashProductLieElem{ZZRingElem}, c::ZZRingElem)
    return parent(e)(data(e) * c)
end

function Base.:*(c::C, e::SmashProductLieElem{C}) where {C <: RingElem}
    coefficient_ring(e) != parent(c) && error("Incompatible rings.")
    return parent(e)(c * data(e))
end

function Base.:*(c::U, e::SmashProductLieElem) where {U <: Union{Rational, IntegerUnion}}
    return parent(e)(c * data(e))
end

function Base.:*(c::ZZRingElem, e::SmashProductLieElem{ZZRingElem})
    return parent(e)(c * data(e))
end

function Base.:^(e::SmashProductLieElem, n::Int)
    return parent(e)(data(e)^n)
end

function comm(e1::SmashProductLieElem, e2::SmashProductLieElem)
    check_parent(e1, e2)
    return parent(e1)(data(e1) * data(e2) - data(e2) * data(e1))
end

function symmetrize(e::SmashProductLieElem)
    return parent(e)(symmetrize(data(e)))
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(
    e1::SmashProductLieElem{C, LieC, LieT},
    e2::SmashProductLieElem{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    return parent(e1) === parent(e2) && data(simplify(e1)) == data(simplify(e2))
end

function Base.hash(e::SmashProductLieElem, h::UInt)
    e = simplify(e)
    b = 0xdcc11ff793ca4ada % UInt
    h = hash(parent(e), h)
    h = hash(data(e), h)
    return xor(h, b)
end

###############################################################################
#
#   Mutable arithmetic
#
###############################################################################

function add!(z::SmashProductLieElem, x::SmashProductLieElem, y::SmashProductLieElem)
    z.alg_elem = add!(z.alg_elem, data(x), data(y))
    z.simplified = false
    return z
end

function sub!(z::SmashProductLieElem, x::SmashProductLieElem, y::SmashProductLieElem)
    z.alg_elem = sub!(z.alg_elem, data(x), data(y))
    z.simplified = false
    return z
end

function mul!(z::SmashProductLieElem{C}, x::SmashProductLieElem{C}, c::Union{Integer, C}) where {C <: RingElem}
    z.alg_elem = mul!(z.alg_elem, data(x), c)
    z.simplified = false
    return z
end

###############################################################################
#
#   Simplification
#
###############################################################################

function simplify(
    e::SmashProductLieElem{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    e.simplified && return e
    e.alg_elem =
        _normal_form(data(e), parent(e).rels::Matrix{Union{Nothing, elem_type(free_associative_algebra_type(C))}})
    e.simplified = true
    return e
end

function _normal_form(a::F, rels::Matrix{Union{Nothing, F}}) where {C <: RingElem, F <: FreeAssociativeAlgebraElem{C}}
    result = zero(parent(a))
    CR = coefficient_ring(a)
    A = parent(a)
    tmp = zero(A)
    while a.length > 0
        c = leading_coefficient(a)
        exp = leading_exponent_word(a)
        t = leading_term(a)
        # 2-arg mutable arithmetic are way slower than 3-arg at the time of writing. TODO: replace once the situation has improved
        tmp = sub!(tmp, a, t)
        a, tmp = tmp, a

        changed = false
        for i in 1:length(exp)-1
            exp[i] > exp[i+1] || continue
            rel = rels[exp[i], exp[i+1]]
            if !isnothing(rel)
                changed = true
                new_term = A([c], [exp[1:i-1]]) * rel * A([one(CR)], [exp[i+2:end]])
                tmp = add!(tmp, a, new_term)
                a, tmp = tmp, a
                break
            end
        end
        if !changed
            tmp = add!(tmp, result, t)
            result, tmp = tmp, result
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
    return smash_product(coefficient_ring(L), L, V)
end

"""
    smash_product(R::Ring, L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElem}

Construct the smash product ``TV \\rtimes U(L)`` and extend the coefficients to `R`.
"""
function smash_product(R::Ring, L::LieAlgebra{C}, V::LieAlgebraModule{C}) where {C <: RingElem}
    @req L == base_lie_algebra(V) "Incompatible module."
    try
        R(zero(coefficient_ring(L)))
    catch e
        if e isa MethodError
            error("Coefficient extension from $(coefficient_ring(L)) to $(R) not possible.")
        else
            rethrow(e)
        end
    end

    dimL = dim(L)
    dimV = dim(V)

    f_alg, _ = free_associative_algebra(
        R,
        [symbols(L); _is_standard_module(V) ? symbols(V) : (x -> Symbol("($x)")).(symbols(V))];
        cached=false,
    )
    f_basisL = [gen(f_alg, i) for i in 1:dimL]
    f_basisV = [gen(f_alg, dimL + i) for i in 1:dimV]

    rels = Matrix{Union{Nothing, elem_type(free_associative_algebra_type(R))}}(nothing, dimL + dimV, dimL + dimV)

    for (i, xi) in enumerate(basis(L)), (j, xj) in enumerate(basis(L))
        commutator =
            sum(R(c) * f_basisL[k] for (k, c) in enumerate(coefficients(xi * xj)) if !iszero(c); init=zero(f_alg))
        rels[i, j] = f_basisL[j] * f_basisL[i] + commutator

    end

    for (i, xi) in enumerate(basis(L)), (j, vj) in enumerate(basis(V))
        commutator =
            sum(R(c) * f_basisV[k] for (k, c) in enumerate(coefficients(xi * vj)) if !iszero(c); init=zero(f_alg))
        rels[i, dimL+j] = f_basisV[j] * f_basisL[i] + commutator
        rels[dimL+j, i] = f_basisL[i] * f_basisV[j] - commutator
    end

    Sp = SmashProductLie{elem_type(R), C, elem_type(L)}(R, L, V, f_alg, rels)

    return Sp
end
