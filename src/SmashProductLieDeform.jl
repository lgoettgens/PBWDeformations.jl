"""
The struct representing a deformation of a Lie algebra smash product.
It consists of the underlying FreeAssAlgebra with relations and some metadata.
It gets created by calling [`deform`](@ref).
"""
@attributes mutable struct SmashProductLieDeform{C <: RingElem, CL <: RingElem} <: NCRing
    sp::SmashProductLie{C, CL}
    rels::QuadraticRelations{C}
    kappa::DeformationMap{C}

    # default constructor for @attributes
    function SmashProductLieDeform{C, CL}(
        sp::SmashProductLie{C, CL},
        rels::QuadraticRelations{C},
        kappa::DeformationMap{C},
    ) where {C <: RingElem, CL <: RingElem}
        new{C, CL}(sp, rels, kappa)
    end
end

mutable struct SmashProductLieDeformElem{C <: RingElem, CL <: RingElem} <: NCRingElem
    p::SmashProductLieDeform{C, CL}   # parent
    alg_elem::FreeAssAlgElem{C}
    simplified::Bool

    function SmashProductLieDeformElem(
        p::SmashProductLieDeform{C, CL},
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

parent_type(::Type{SmashProductLieDeformElem{C, CL}}) where {C <: RingElem, CL <: RingElem} =
    SmashProductLieDeform{C, CL}

elem_type(::Type{SmashProductLieDeform{C, CL}}) where {C <: RingElem, CL <: RingElem} = SmashProductLieDeformElem{C, CL}

parent(e::SmashProductLieDeformElem) = e.p

base_ring(D::SmashProductLieDeform) = base_ring(lie_algebra(D))

base_ring(e::SmashProductLieDeformElem) = base_ring(parent(e))

lie_algebra(D::SmashProductLieDeform) = lie_algebra(D.sp)

lie_module(D::SmashProductLieDeform) = lie_module(D.sp)

underlying_algebra(D::SmashProductLieDeform) = underlying_algebra(D.sp) # TODO: create new algebra for D

ngens(D::SmashProductLieDeform) = ngens(underlying_algebra(D))
function ngens(D::SmashProductLieDeform, part::Symbol)
    part == :L && return dim(lie_algebra(D))
    part == :V && return dim(lie_module(D))
    error("Invalid part.")
end

gens(D::SmashProductLieDeform) = map(D, gens(underlying_algebra(D)))
function gens(D::SmashProductLieDeform, part::Symbol)
    part == :L && return [D(gen(underlying_algebra(D), i)) for i in 1:dim(lie_algebra(D))]
    part == :V && return [D(gen(underlying_algebra(D), i + dim(lie_algebra(D)))) for i in 1:dim(lie_module(D))]
    error("Invalid part.")
end

gen(D::SmashProductLieDeform, i::Int) = D(gen(underlying_algebra(D), i))
function gen(D::SmashProductLieDeform, i::Int, part::Symbol)
    @req 1 <= i <= ngens(D, part) "Invalid generator index."
    part == :L && return D(gen(underlying_algebra(D), i))
    part == :V && return D(gen(underlying_algebra(D), i + dim(lie_algebra(D))))
    error("Invalid part.")
end

function zero(D::SmashProductLieDeform)
    return D(zero(underlying_algebra(D)))
end

function iszero(e::SmashProductLieDeformElem)
    return iszero(simplify!(e).alg_elem)
end

function one(D::SmashProductLieDeform)
    return D(one(underlying_algebra(D)))
end

function isone(e::SmashProductLieDeformElem)
    return isone(simplify!(e).alg_elem)
end

function Base.deepcopy_internal(e::SmashProductLieDeformElem, dict::IdDict)
    return SmashProductLieDeformElem(parent(e), deepcopy_internal(e.alg_elem, dict); simplified=e.simplified)
end

function check_parent(e1::SmashProductLieDeformElem{C}, e2::SmashProductLieDeformElem{C}) where {C <: RingElem}
    parent(e1) != parent(e2) && error("Incompatible smash product deformations.")
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

function (D::SmashProductLieDeform{C, CL})(e::FreeAssAlgElem{C}) where {C <: RingElem, CL <: RingElem}
    if underlying_algebra(D) !== parent(e)
        e = underlying_algebra(D)(e)
    end
    return SmashProductLieDeformElem(D, e)
end

function (D::SmashProductLieDeform{C, CL})(e::SmashProductLieDeformElem{C, CL}) where {C <: RingElem, CL <: RingElem}
    @req parent(e) == D "Incompatible smash product deformations."
    return e
end

function (D::SmashProductLieDeform{C, CL})(e::SmashProductLieElem{C, CL}) where {C <: RingElem, CL <: RingElem}
    @req parent(e) == D.sp "Incompatible smash products."
    return D(0) # TODO
end

###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(e::SmashProductLieDeformElem)
    return parent(e)(-e.alg_elem)
end

function Base.:+(e1::SmashProductLieDeformElem{C}, e2::SmashProductLieDeformElem{C}) where {C <: RingElem}
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem + e2.alg_elem)
end

function Base.:-(e1::SmashProductLieDeformElem{C}, e2::SmashProductLieDeformElem{C}) where {C <: RingElem}
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem - e2.alg_elem)
end

function Base.:*(e1::SmashProductLieDeformElem{C}, e2::SmashProductLieDeformElem{C}) where {C <: RingElem}
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem * e2.alg_elem)
end

function Base.:*(e::SmashProductLieDeformElem{C}, c::C) where {C <: RingElem}
    base_ring(e) != parent(c) && error("Incompatible rings.")
    return parent(e)(e.alg_elem * c)
end

function Base.:*(e::SmashProductLieDeformElem{C}, c::U) where {C <: RingElem, U <: Union{Rational, Integer}}
    return parent(e)(e.alg_elem * c)
end

function Base.:*(c::C, e::SmashProductLieDeformElem{C}) where {C <: RingElem}
    base_ring(e) != parent(c) && error("Incompatible rings.")
    return parent(e)(c * e.alg_elem)
end

function Base.:*(c::U, e::SmashProductLieDeformElem{C}) where {C <: RingElem, U <: Union{Rational, Integer}}
    return parent(e)(c * e.alg_elem)
end

function Base.:^(e::SmashProductLieDeformElem, n::Int)
    return parent(e)(e.alg_elem^n)
end

function comm(e1::SmashProductLieDeformElem{C}, e2::SmashProductLieDeformElem{C}) where {C <: RingElem}
    check_parent(e1, e2)
    return parent(e1)(e1.alg_elem * e2.alg_elem - e2.alg_elem * e1.alg_elem)
end

###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(e1::SmashProductLieDeformElem, e2::SmashProductLieDeformElem)
    return parent(e1) === parent(e2) && simplify!(e1).alg_elem == simplify!(e2).alg_elem
end

function Base.hash(e::SmashProductLieDeformElem, h::UInt)
    e = simplify!(e)
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

function simplify!(e::SmashProductLieDeformElem)
    e.simplified && return e
    e.alg_elem = normal_form(e.alg_elem, parent(e).rels)
    e.simplified = true
    return e
end

function simplify(e::SmashProductLieDeformElem)
    return deepcopy(e) |> simplify!
end

###############################################################################
#
#   Constructor
#
###############################################################################

"""
    deform(sp::SmashProductLie{C}, kappa::DeformationMap{C}) where {C <: RingElem}

Constructs the deformation of the smash product `sp` by the deformation map `kappa`.

Returns a [`SmashProductDeformLie`](@ref) struct and a two-part basis.
"""
function deform(sp::SmashProductLie{C, CL}, kappa::DeformationMap{C}) where {C <: RingElem, CL <: RingElem}
    dimL = dim(sp.L)
    dimV = dim(sp.V)

    @req size(kappa) == (dimV, dimV) "kappa has wrong dimensions."

    basisV = [gen(sp.alg, dimL + i) for i in 1:dimV]

    for i in 1:dimV, j in 1:i
        @req kappa[i, j] == -kappa[j, i] "kappa is not skew-symmetric."
        @req all(x -> x <= dimL, Iterators.flatten(exponent_words(kappa[i, j]))) "kappa does not only take values in the hopf algebra"
        @req all(x -> x <= dimL, Iterators.flatten(exponent_words(kappa[j, i]))) "kappa does not only take values in the hopf algebra"
    end

    symmetric = true
    rels = deepcopy(sp.rels)
    for i in 1:dimV, j in 1:dimV
        # We have the commutator relation [v_i, v_j] = kappa[i,j]
        # which is equivalent to v_i*v_j = v_j*v_i + kappa[i,j]
        rels[(dimL + i, dimL + j)] = basisV[j] * basisV[i] + kappa[i, j]
        symmetric &= iszero(kappa[i, j])
    end

    d = SmashProductLieDeform{C, CL}(sp, rels, kappa)

    set_attribute!(d, :is_symmetric, symmetric)

    return d
end

"""
    symmetric_deformation(sp::SmashProductLie{C}) where {C <: RingElem}

Constructs the symmetric deformation of the smash product `sp`.
"""
function symmetric_deformation(sp::SmashProductLie{C, CL}) where {C <: RingElem, CL <: RingElem}
    kappa = fill(zero(sp.alg), dim(sp.V), dim(sp.V))
    d = deform(sp, kappa)
    return d
end
