"""
Concrete subtype of [`DeformBasis`](@ref) that implements the standard basis.
Each element of the basis is a skew-symmetric matrix with 2 non-zero entries,
where one entry is a pure tensor power of degree ∈ `degs` over the Lie algebra part
of the smash product, and the other entry is its additive inverse.
"""
struct StdDeformBasis{T <: SmashProductLieElem} <: DeformBasis{T}
    len::Int
    iter
    extra_data::Dict{DeformationMap{T}, Nothing}
    normalize

    function StdDeformBasis{T}(
        sp::SmashProductLie{C, LieC, LieT},
        degs::AbstractVector{Int}
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}, T <: SmashProductLieElem{C, LieC, LieT}}
        dimL = dim(base_lie_algebra(sp))
        dimV = dim(base_module(sp))
        iter = (
            begin
                kappa = zero_matrix(sp, dimV, dimV)
                entry = prod(map(k -> gen(sp, k, :L), ind); init=one(sp))
                kappa[i, j] += entry
                kappa[j, i] -= entry
                kappa
            end for i in 1:dimV for j in i+1:dimV for d in degs for ind in multicombinations(1:dimL, d)
        )

        len = div(dimV * (dimV - 1), 2) * sum(binomial(dimL + k - 1, k) for k in degs)
        return new{T}(len, iter, Dict{DeformationMap{T}, Nothing}(), normalize_default)
    end
end

function Base.iterate(i::StdDeformBasis)
    return iterate(i.iter)
end

function Base.iterate(i::StdDeformBasis, s)
    return iterate(i.iter, s)
end

Base.length(basis::StdDeformBasis) = basis.len
