"""
Concrete subtype of [`DeformBasis`](@ref) that implements the standard basis.
Each element of the basis is a skew-symmetric matrix with 2 non-zero entries,
where one entry is a pure tensor power of degree ∈ `degs` over the Lie algebra part
of the smash product, and the other entry is its additive inverse.
"""
struct StdDeformBasis{T <: SmashProductLieElem} <: DeformBasis{T}
    sp::SmashProductLie # parent_type(T)
    degs::Vector{Int}
    len::Int
    iter

    function StdDeformBasis(sp::SmashProductLie, degs::AbstractVector{Int})
        return StdDeformBasis(sp, collect(degs))
    end

    function StdDeformBasis(sp::SmashProductLie, degs::Vector{Int})
        @req coefficient_ring(sp) === coefficient_ring(base_lie_algebra(sp)) "Deformation bases don't support extension of the coefficient ring of the smash product."
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
        return StdDeformBasis{elem_type(sp)}(sp, degs, len, iter)
    end

    function StdDeformBasis{T}(sp::SmashProductLie, degs::Vector{Int}, len::Int, iter) where {T <: SmashProductLieElem}
        @assert sp isa parent_type(T)
        # This inner constructor just sets the fields directly, without any checks.
        return new{T}(sp, degs, len, iter)
    end
end

function Base.iterate(i::StdDeformBasis)
    return iterate(i.iter)::Union{Tuple{eltype(i), Any}, Nothing}
end

function Base.iterate(i::StdDeformBasis, s)
    return iterate(i.iter, s)::Union{Tuple{eltype(i), Any}, Nothing}
end

Base.length(basis::StdDeformBasis) = basis.len

function Base.show(io::IO, basis::StdDeformBasis)
    print(terse(pretty(io)), "StdDeformBasis of ", Lowercase(), basis.sp, " with degrees ", basis.degs)
end
