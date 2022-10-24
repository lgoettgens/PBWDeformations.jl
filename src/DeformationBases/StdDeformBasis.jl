"""
Concrete subtype of [`DeformBasis`](@ref) that implements the standard basis.
Each element of the basis is a skew-symmetric matrix with 2 non-zero entries,
where one entry is a pure tensor power of degree ∈ `degs` over the Lie algebra part
of the smash product, and the other entry is its additive inverse.
"""
struct StdDeformBasis{C <: RingElement} <: DeformBasis{C}
    len::Int
    iter

    function StdDeformBasis{C}(sp::SmashProductLie{C}, degs::AbstractVector{Int}) where {C <: RingElement}
        dimL = sp.dimL
        dimV = sp.dimV
        R = coefficient_ring(sp.alg)
        iter = (
            begin
                kappa = fill(sp.alg(0), dimV, dimV)
                entry = prod(map(k -> sp.basisL[k], ind); init=sp.alg(1))
                kappa[i, j] += entry
                kappa[j, i] -= entry
                kappa
            end for i in 1:dimV for j in i+1:dimV for d in degs for
            ind in Combinatorics.with_replacement_combinations(1:dimL, d)
        )

        len = div(dimV * (dimV - 1), 2) * sum(binomial(dimL + k - 1, k) for k in degs)
        return new{C}(len, iter)
    end
end

function Base.iterate(i::StdDeformBasis)
    return iterate(i.iter)
end

function Base.iterate(i::StdDeformBasis, s)
    return iterate(i.iter, s)
end

Base.length(base::StdDeformBasis) = base.len
