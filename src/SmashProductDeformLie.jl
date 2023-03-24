"""
The struct representing a deformation of a Lie algebra smash product.
It consists of the underlying FreeAssAlgebra with relations and some metadata.
It gets created by calling [`deform`](@ref).
"""
@attributes mutable struct SmashProductDeformLie{C <: RingElement, C_lie <: RingElement}
    sp::SmashProductLie{C, C_lie}
    rels::QuadraticRelations{C}
    kappa::DeformationMap{C}

    # default constructor for @attributes
    function SmashProductDeformLie{C, C_lie}(
        sp::SmashProductLie{C, C_lie},
        rels::QuadraticRelations{C},
        kappa::DeformationMap{C},
    ) where {C <: RingElement, C_lie <: RingElement}
        new{C, C_lie}(sp, rels, kappa)
    end
end


"""
    deform(sp::SmashProductLie{C}, kappa::DeformationMap{C}) where {C <: RingElement}

Constructs the deformation of the smash product `sp` by the deformation map `kappa`.

Returns a [`SmashProductDeformLie`](@ref) struct and a two-part basis.
"""
function deform(sp::SmashProductLie{C, C_lie}, kappa::DeformationMap{C}) where {C <: RingElement, C_lie <: RingElement}
    dimL = dim(sp.L)
    dimV = dim(sp.V)

    size(kappa) == (dimV, dimV) || throw(ArgumentError("kappa has wrong dimensions."))

    basisV = [gen(sp.alg, dimL + i) for i in 1:dimV]

    for i in 1:dimV, j in 1:i
        kappa[i, j] == -kappa[j, i] || throw(ArgumentError("kappa is not skew-symmetric."))
        all(x -> x <= dimL, Iterators.flatten(exponent_words(kappa[i, j]))) ||
            throw(ArgumentError("kappa does not only take values in the hopf algebra"))
        all(x -> x <= dimL, Iterators.flatten(exponent_words(kappa[j, i]))) ||
            throw(ArgumentError("kappa does not only take values in the hopf algebra"))
    end

    symmetric = true
    rels = deepcopy(sp.rels)
    for i in 1:dimV, j in 1:dimV
        # We have the commutator relation [v_i, v_j] = kappa[i,j]
        # which is equivalent to v_i*v_j = v_j*v_i + kappa[i,j]
        rels[(dimL + i, dimL + j)] = basisV[j] * basisV[i] + kappa[i, j]
        symmetric &= iszero(kappa[i, j])
    end

    d = SmashProductDeformLie{C, C_lie}(sp, rels, kappa)

    set_attribute!(d, :is_symmetric, symmetric)

    return d
end

"""
    symmetric_deformation(sp::SmashProductLie{C}) where {C <: RingElement}

Constructs the symmetric deformation of the smash product `sp`.
"""
function symmetric_deformation(sp::SmashProductLie{C, C_lie}) where {C <: RingElement, C_lie <: RingElement}
    kappa = fill(zero(sp.alg), dim(sp.V), dim(sp.V))
    d = deform(sp, kappa)
    return d
end


ngens(d::SmashProductDeformLie{C}) where {C <: RingElement} = length(gens(d.sp.alg)) # ngens(d.sp.alg), see https://github.com/Nemocas/AbstractAlgebra.jl/pull/1295
function ngens(d::SmashProductDeformLie{C}, part::Symbol) where {C <: RingElement}
    part == :L && return dim(d.sp.L)
    part == :V && return dim(d.sp.V)
    error("Invalid part.")
end

gens(d::SmashProductDeformLie{C}) where {C <: RingElement} = gens(d.sp.alg)
function gens(d::SmashProductDeformLie{C}, part::Symbol) where {C}
    part == :L && return [gen(d.sp.alg, i) for i in 1:dim(d.sp.L)]
    part == :V && return [gen(d.sp.alg, i + dim(d.sp.L)) for i in 1:dim(d.sp.V)]
    error("Invalid part.")
end

gen(d::SmashProductDeformLie{C}, i::Int) where {C <: RingElement} = gen(d.sp.alg, i)
function gen(d::SmashProductDeformLie{C}, i::Int, part::Symbol) where {C <: RingElement}
    1 <= i <= ngens(d, part) || error("Invalid generator index.")
    part == :L && return gen(d.sp.alg, i)
    part == :V && return gen(d.sp.alg, i + dim(d.sp.L))
    error("Invalid part.")
end


function show(io::IO, d::SmashProductDeformLie)
    if get_attribute(d, :is_symmetric, false)
        print(io, "Symmetric deformation of ")
    else
        print(io, "Deformation of ")
    end
    print(IOContext(io, :compact => true), d.sp)
end

function change_base_ring(R::Ring, d::SmashProductDeformLie{C}) where {C <: RingElement}
    sp = change_base_ring(R, d.sp)
    rels = QuadraticRelations{elem_type(R)}(k => change_base_ring(R, a, parent=alg) for (k, a) in d.rels)
    kappa = map(e -> change_base_ring(R, e, parent=alg), d.kappa)

    return SmashProductDeformLie{elem_type(R)}(sp, rels, kappa)
end
