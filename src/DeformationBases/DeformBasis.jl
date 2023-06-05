"""
    DeformationMap{C} = Matrix{FreeAssAlgElem{C}} where {C <: RingElem}

The type for deformation maps of a Lie algebra smash product.
The entry `kappa[i,j]` should be the image of ``v_i \\wedge v_j`` under the deformation map, i.e. ``κ(v_i,v_j)``.
Deformation maps are always assumed to be quadratic and skew-symmetric.
"""
DeformationMap{C} = Matrix{<:FreeAssAlgElem{C}} where {C <: RingElem}


"""
    abstract type DeformBasis{C <: RingElem} end

A basis for a deformation map space of a Lie algebra smash product.
The constructor of a subtype should accept a [`SmashProductLie`](@ref) and an `AbstractVector{Int}` of degrees.
It is required that `Base.length` and `Base.iterate` are implemented for subtypes,
where iterating yields objects of type `DeformationMap{C}`.

For a reference implementation, we refer to [`StdDeformBasis`](@ref).
"""
abstract type DeformBasis{C <: RingElem} end

Base.eltype(::Type{DeformBasis{C}}) where {C <: RingElem} = DeformationMap{C}

Base.length(basis::DeformBasis) = error("length not implemented for $(typeof(basis))")


"""
    lookup_data(m::DeformationMap{C}, basis::DeformBasis{C}) where {C <: RingElem}

Look up additional data that was used to generate the deformation map `m` in the basis `basis`.
This can e.g. be an arc diagram or a pseudograph.
"""
function lookup_data(m::DeformationMap{C}, basis::DeformBasis{C}) where {C <: RingElem}
    m = basis.normalize(m)
    if haskey(basis.extra_data, m)
        return basis.extra_data[m]
    else
        return nothing
    end
end

function normalize_default(m::DeformationMap{C}) where {C <: RingElem}
    nz_index = findfirst(x -> !iszero(x), m)
    if nz_index === nothing
        return m
    end
    cu = canonical_unit(m[CartesianIndex(nz_index[2], nz_index[1])])
    m = map(e -> C(1 // cu) * e, m)
end
