Base.eltype(::Type{<:DeformBasis{T}}) where {T <: SmashProductLieElem} = DeformationMap{T}

Base.length(basis::DeformBasis) = error("length not implemented for $(typeof(basis))")

function normalize_scaling(m::DeformationMap{T}) where {T <: SmashProductLieElem}
    nz_index = findfirst(x -> !iszero(x), m) # the index is in the lower triangular part
    if isnothing(nz_index)
        return m, zero(coefficient_ring(base_ring(m)))
    end
    coeff = -leading_coefficient(data(m[nz_index])) # negative to change to the upper triangular part
    coeff_inv = inv(coeff)
    m = map_entries(e -> coeff_inv * e, m)
    return m, coeff
end

function filter_independent(R, input)
    deduplicated = unique(Iterators.filter(!iszero, input))::Vector{eltype(input)}
    _, rels = is_linearly_independent_with_relations(R, deduplicated)
    inds = [findlast(!iszero, vec(rels[i, :]))::Int for i in 1:nrows(rels)]
    deleteat!(deduplicated, inds)
    return deduplicated, length(deduplicated)
end
