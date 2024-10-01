###############################################################################
#
#   Basic manipulation
#
###############################################################################

arent_type(
    ::Type{SmashProductLieDeformElem{C, LieC, LieT}},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = SmashProductLieDeform{C, LieC, LieT}

elem_type(
    ::Type{SmashProductLieDeform{C, LieC, LieT}},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = SmashProductLieDeformElem{C, LieC, LieT}

parent(e::SmashProductLieDeformElem) = e.p

coefficient_ring(
    D::SmashProductLieDeform{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = coefficient_ring(D.sp)::parent_type(C)

coefficient_ring(e::SmashProductLieDeformElem) = coefficient_ring(parent(e))

base_lie_algebra(
    D::SmashProductLieDeform{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = base_lie_algebra(D.sp)::parent_type(LieT)

base_module(
    D::SmashProductLieDeform{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}} = base_module(D.sp)::LieAlgebraModule{LieC}

underlying_algebra(D::SmashProductLieDeform) = underlying_algebra(D.sp) # TODO: create new algebra for D

ngens(D::SmashProductLieDeform) = ngens(underlying_algebra(D))
function ngens(D::SmashProductLieDeform, part::Symbol)
    part == :L && return dim(base_lie_algebra(D))
    part == :V && return dim(base_module(D))
    error("Invalid part.")
end

gens(D::SmashProductLieDeform) = map(D, gens(underlying_algebra(D)))
function gens(D::SmashProductLieDeform, part::Symbol)
    part == :L && return [D(gen(underlying_algebra(D), i)) for i in 1:dim(base_lie_algebra(D))]
    part == :V && return [D(gen(underlying_algebra(D), i + dim(base_lie_algebra(D)))) for i in 1:dim(base_module(D))]
    error("Invalid part.")
end

gen(D::SmashProductLieDeform, i::Int) = D(gen(underlying_algebra(D), i))
function gen(D::SmashProductLieDeform, i::Int, part::Symbol)
    @req 1 <= i <= ngens(D, part) "Invalid generator index."
    part == :L && return D(gen(underlying_algebra(D), i))
    part == :V && return D(gen(underlying_algebra(D), i + dim(base_lie_algebra(D))))
    error("Invalid part.")
end

function zero(D::SmashProductLieDeform)
    return D(zero(underlying_algebra(D)))
end

function iszero(e::SmashProductLieDeformElem)
    return iszero(simplify(e).alg_elem)
end

function one(D::SmashProductLieDeform)
    return D(one(underlying_algebra(D)))
end

function isone(e::SmashProductLieDeformElem)
    return isone(simplify(e).alg_elem)
end

function Base.deepcopy_internal(e::SmashProductLieDeformElem, dict::IdDict)
    return SmashProductLieDeformElem(parent(e), deepcopy_internal(e.alg_elem, dict); simplified=e.simplified)
end

function check_parent(
    e1::SmashProductLieDeformElem{C, LieC, LieT},
    e2::SmashProductLieDeformElem{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    parent(e1) != parent(e2) && error("Incompatible smash product deformations.")
end

function AbstractAlgebra.promote_rule(
    T1::Type{<:SmashProductLieDeformElem{C1}},
    ::Type{C2},
) where {C1 <: RingElem, C2 <: RingElem}
    AbstractAlgebra.promote_rule(C1, C2) == C1 ? T1 : Union{}
end

###############################################################################
#
#   String I/O
#
###############################################################################

function show(io::IO, D::SmashProductLieDeform)
    if get_attribute(D, :is_symmetric, false)
        print(io, "Symmetric deformation of ")
    else
        print(io, "Deformation of ")
    end
    print(IOContext(io, :compact => true), D.sp)
end


function show(io::IO, e::SmashProductLieDeformElem)
    show(io, e.alg_elem)
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (D::SmashProductLieDeform)()
    return zero(D)
end

function (D::SmashProductLieDeform)(e::Union{RingElement, NCRingElem})
    return D(underlying_algebra(D)(e))
end

function (D::SmashProductLieDeform{C})(e::FreeAssAlgElem{C}) where {C <: RingElem}
    if underlying_algebra(D) !== parent(e)
        e = underlying_algebra(D)(e)
    end
    return SmashProductLieDeformElem(D, e)
end

function (D::SmashProductLieDeform{C, LieC, LieT})(
    e::SmashProductLieDeformElem{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req parent(e) == D "Incompatible smash product deformations."
    return e
end

function (D::SmashProductLieDeform{C, LieC, LieT})(
    e::SmashProductLieElem{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req parent(e) == D.sp "Incompatible smash products."
    return D(e.alg_elem)
end

function (D::SmashProductLieDeform{C, LieC, LieT})(
    x::LieT,
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req parent(x) == base_lie_algebra(D) "Incompatible smash products."
    return sum(c * b for (c, b) in zip(coefficients(x), gens(D, :L)))
end

function (D::SmashProductLieDeform{C, LieC, LieT})(
    v::LieAlgebraModuleElem{LieC},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req parent(v) == base_module(D) "Incompatible smash products."
    return sum(c * b for (c, b) in zip(coefficients(v), gens(D, :V)))
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(e::SmashProductLieDeformElem)
    return parent(e)(-e.alg_elem)
end

function Base.:+(e1::SmashProductLieDeformElem, e2::SmashProductLieDeformElem)
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem + e2.alg_elem)
end

function Base.:-(e1::SmashProductLieDeformElem, e2::SmashProductLieDeformElem)
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem - e2.alg_elem)
end

function Base.:*(e1::SmashProductLieDeformElem, e2::SmashProductLieDeformElem)
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem * e2.alg_elem)
end

function Base.:*(e::SmashProductLieDeformElem{C}, c::C) where {C <: RingElem}
    coefficient_ring(e) != parent(c) && error("Incompatible rings.")
    return parent(e)(e.alg_elem * c)
end

function Base.:*(e::SmashProductLieDeformElem, c::U) where {U <: Union{Rational, IntegerUnion}}
    return parent(e)(e.alg_elem * c)
end

function Base.:*(e::SmashProductLieDeformElem{ZZRingElem}, c::ZZRingElem)
    return parent(e)(e.alg_elem * c)
end

function Base.:*(c::C, e::SmashProductLieDeformElem{C}) where {C <: RingElem}
    coefficient_ring(e) != parent(c) && error("Incompatible rings.")
    return parent(e)(c * e.alg_elem)
end

function Base.:*(c::U, e::SmashProductLieDeformElem) where {U <: Union{Rational, IntegerUnion}}
    return parent(e)(c * e.alg_elem)
end

function Base.:*(c::ZZRingElem, e::SmashProductLieDeformElem{ZZRingElem})
    return parent(e)(c * e.alg_elem)
end

function Base.:^(e::SmashProductLieDeformElem, n::Int)
    return parent(e)(e.alg_elem^n)
end

function comm(e1::SmashProductLieDeformElem, e2::SmashProductLieDeformElem)
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem * e2.alg_elem - e2.alg_elem * e1.alg_elem)
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(
    e1::SmashProductLieDeformElem{C, LieC, LieT},
    e2::SmashProductLieDeformElem{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    return parent(e1) === parent(e2) && simplify(e1).alg_elem == simplify(e2).alg_elem
end

function Base.hash(e::SmashProductLieDeformElem, h::UInt)
    e = simplify(e)
    b = 0x97eb07aa70e4a59c % UInt
    h = hash(parent(e), h)
    h = hash(e.alg_elem, h)
    return xor(h, b)
end

###############################################################################
#
#   Simplification
#
###############################################################################

function simplify(e::SmashProductLieDeformElem)
    e.simplified && return e
    e.alg_elem = _normal_form(e.alg_elem, parent(e).rels)
    e.simplified = true
    return e
end

###############################################################################
#
#   Constructor
#
###############################################################################

"""
    deform(sp::SmashProductLie{C}, kappa::DeformationMap{elem_type(sp)}) where {C <: RingElem}

Constructs the deformation of the smash product `sp` by the deformation map `kappa`.

Returns a [`SmashProductLieDeform`](@ref) struct and a two-part basis.
"""
function deform(
    sp::SmashProductLie{C, LieC, LieT},
    kappa::DeformationMap{SmashProductLieElem{C, LieC, LieT}},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    dimL = dim(base_lie_algebra(sp))
    dimV = dim(base_module(sp))

    @req size(kappa) == (dimV, dimV) "kappa has wrong dimensions."
    @req all(e -> parent(e) == sp, kappa) "Incompatible smash products."

    basisV = [gen(underlying_algebra(sp), dimL + i) for i in 1:dimV]

    for i in 1:dimV, j in 1:i
        @req kappa[i, j] == -kappa[j, i] "kappa is not skew-symmetric."
        @req all(<=(dimL), Iterators.flatten(exponent_words(kappa[i, j].alg_elem))) "kappa does not only take values in the hopf algebra"
        @req all(<=(dimL), Iterators.flatten(exponent_words(kappa[j, i].alg_elem))) "kappa does not only take values in the hopf algebra"
    end

    symmetric = true
    rels = deepcopy(sp.rels)
    for i in 1:dimV, j in 1:dimV
        # We have the commutator relation [v_i, v_j] = kappa[i,j]
        # which is equivalent to v_i*v_j = v_j*v_i + kappa[i,j]
        rels[dimL+i, dimL+j] = basisV[j] * basisV[i] + simplify(kappa[i, j]).alg_elem
        symmetric &= iszero(kappa[i, j])
    end

    d = SmashProductLieDeform{C, LieC, LieT}(sp, rels, kappa)

    set_attribute!(d, :is_symmetric, symmetric)

    return d
end

function deform(
    sp::SmashProductLie{C, LieC, LieT},
    kappa::MatElem{<:FreeAssAlgElem{C}},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    return deform(sp, map_entries(sp, kappa))
end

"""
    symmetric_deformation(sp::SmashProductLie{C}) where {C <: RingElem}

Constructs the symmetric deformation of the smash product `sp`.
"""
function symmetric_deformation(
    sp::SmashProductLie{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    kappa = zero_matrix(sp, dim(base_module(sp)), dim(base_module(sp)))
    d = deform(sp, kappa)
    return d
end
