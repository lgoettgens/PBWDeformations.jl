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


function normalize_basis(basiselems)
    unique((
        begin
            first_nz = begin
                ind = findfirst(x -> !iszero(x), b)
                CartesianIndex(ind[2], ind[1])
            end
            cu = canonical_unit(b[first_nz])
            b = map(e -> divexact(e, cu), b)
        end
    ) for b in basiselems if !iszero(b))
end
