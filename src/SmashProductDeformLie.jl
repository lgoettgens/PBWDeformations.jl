"""
The struct representing a deformation of a Lie algebra smash product.
It consists of the underlying FreeAssAlgebra with relations and some metadata.
It gets created by calling [`smash_product_deform_lie`](@ref).
"""
mutable struct SmashProductDeformLie{C <: RingElement}
    dimL::Int
    dimV::Int
    basisL::Vector{FreeAssAlgElem{C}}
    basisV::Vector{FreeAssAlgElem{C}}
    coeff_ring::Ring
    alg::FreeAssAlgebra{C}
    rels::QuadraticRelations{C}
    symmetric::Bool
    kappa::DeformationMap{C}
end


"""
    smash_product_deform_lie(sp::SmashProductLie{C}, kappa::DeformationMap{C}) where {C <: RingElement}

Constructs the deformation of the smash product `sp` by the deformation map `kappa`.

Returns a [`SmashProductDeformLie`](@ref) struct and a two-part basis.
"""
function smash_product_deform_lie(sp::SmashProductLie{C}, kappa::DeformationMap{C}) where {C <: RingElement}
    size(kappa) == (sp.dimV, sp.dimV) || throw(ArgumentError("kappa has wrong dimensions."))

    dimL = sp.dimL
    dimV = sp.dimV
    coeff_ring = sp.coeff_ring
    basisL = sp.basisL
    basisV = sp.basisV

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

    return SmashProductDeformLie{C}(dimL, dimV, basisL, basisV, coeff_ring, sp.alg, rels, symmetric, kappa),
    (basisL, basisV)
end

"""
    smash_product_symmdeform_lie(sp::SmashProductLie{C}) where {C <: RingElement}

Constructs the symmetric deformation of the smash product `sp`.
"""
function smash_product_symmdeform_lie(sp::SmashProductLie{C}) where {C <: RingElement}
    kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
    return smash_product_deform_lie(sp, kappa)
end

ngens(d::SmashProductDeformLie) = d.dimL, d.dimV

function gens(d::SmashProductDeformLie{C}) where {C <: RingElement}
    return [gen(d.alg, i) for i in 1:d.dimL], [gen(d.alg, i + d.dimL) for i in 1:d.dimV]
end


function show(io::IO, deform::SmashProductDeformLie)
    local max_gens = 4 # largest number of generators to print
    if deform.symmetric
        print(io, "Symmetric ")
    end
    print(io, "Deformation of ")
    print(io, "Lie Algebra Smash Product with basis ")
    for i in 1:min(deform.dimL - 1, max_gens - 1)
        print(io, string(deform.alg.S[i]), ", ")
    end
    if deform.dimL > max_gens
        print(io, "..., ")
    end
    print(io, string(deform.alg.S[deform.dimL]) * ", ")
    for i in 1:min(deform.dimV - 1, max_gens - 1)
        print(io, string(deform.alg.S[deform.dimL+i]), ", ")
    end
    if deform.dimV > max_gens
        print(io, "..., ")
    end
    print(io, string(deform.alg.S[deform.dimL+deform.dimV]))
    print(io, " over ")
    print(IOContext(io, :compact => true), deform.coeff_ring)
end

function change_base_ring(R::Ring, d::SmashProductDeformLie{C}) where {C <: RingElement}
    alg, _ = FreeAssociativeAlgebra(R, d.alg.S)
    basisL = map(b -> change_base_ring(R, b, parent=alg), d.basisL)
    basisV = map(b -> change_base_ring(R, b, parent=alg), d.basisV)
    kappa = map(e -> change_base_ring(R, e, parent=alg), d.kappa)
    rels = QuadraticRelations{elem_type(R)}(k => change_base_ring(R, a, parent=alg) for (k, a) in d.rels)

    return SmashProductDeformLie{elem_type(R)}(d.dimL, d.dimV, basisL, basisV, R, alg, rels, d.symmetric, kappa)
end
