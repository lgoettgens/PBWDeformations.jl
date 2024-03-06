"""
Concrete subtype of [`DeformBasis`](@ref) that implements the union of multiple deformation bases.
"""
struct UnionBasis{C <: RingElem} <: DeformBasis{C}
    len::Int
    iter
    extra_data::Dict{DeformationMap{C}, Nothing}
    normalize

    function UnionBasis{C}(bases::Vector{<:DeformBasis{C}}) where {C <: RingElem}
        len = sum(length(b) for b in bases; init=0)
        iter = Iterators.flatten(b.iter for b in bases)
        return new{C}(len, iter, Dict{DeformationMap{C}, Nothing}(), normalize_default)
    end
end

function Base.iterate(i::UnionBasis)
    return iterate(i.iter)
end

function Base.iterate(i::UnionBasis, s)
    return iterate(i.iter, s)
end

Base.length(basis::UnionBasis) = basis.len
