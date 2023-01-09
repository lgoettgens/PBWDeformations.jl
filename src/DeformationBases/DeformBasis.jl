"""
    DeformationMap{C} = Matrix{QuadraticQuoAlgebraElem{C}} where {C <: RingElement}

The type for deformation maps of a Lie algebra smash product.
The entry `kappa[i,j]` should be the image of ``v_i \\wedge v_j`` under the deformation map, i.e. ``κ(v_i,v_j)``.
Deformation maps are always assumed to be quadratic and skew-symmetric.
"""
DeformationMap{C} = Matrix{QuadraticQuoAlgebraElem{C}} where {C <: RingElement}


"""
    abstract type DeformBasis{C <: RingElement} end

A basis for a deformation map space of a Lie algebra smash product.
The constructor of a subtype should accept a [`SmashProductLie`](@ref) and an `AbstractVector{Int}` of degrees.
It is required that `Base.length` and `Base.iterate` are implemented for subtypes,
where iterating yields objects of type `DeformationMap{C}`.

For a reference implementation, we refer to [`StdDeformBasis`](@ref).
"""
abstract type DeformBasis{C <: RingElement} end

Base.eltype(::Type{DeformBasis{C}}) where {C <: RingElement} = DeformationMap{C}

function Base.length(base::DeformBasis)
    error("length not implemented for $(typeof(base))")
end

function lookup_data(m::DeformationMap{C}, base::DeformBasis{C}) where {C <: RingElement}
    m = base.normalize(m)
    if haskey(base.extra_data, m)
        return base.extra_data[m]
    else
        return nothing
    end
end

function normalize_default(m::DeformationMap{C}) where {C <: RingElement}
    nz_index = findfirst(x -> !iszero(x), m)
    if nz_index === nothing
        return m
    end
    cu = canonical_unit(m[CartesianIndex(nz_index[2], nz_index[1])])
    m = map(e -> divexact(e, cu), m)
end
