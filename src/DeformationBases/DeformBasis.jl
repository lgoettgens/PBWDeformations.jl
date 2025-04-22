Base.eltype(::Type{<:DeformBasis{T}}) where {T <: SmashProductLieElem} = DeformationMap{T}

Base.length(basis::DeformBasis) = error("length not implemented for $(typeof(basis))")


"""
    lookup_data(m::DeformationMap{T}, basis::DeformBasis{T}) where {T <: SmashProductLieElem}

Look up additional data that was used to generate the deformation map `m` in the basis `basis`.
This can e.g. be an arc diagram or a pseudograph.
"""
function lookup_data(m::DeformationMap{T}, basis::DeformBasis{T}) where {T <: SmashProductLieElem}
    if !basis.no_normalize
        m = normalize(m)
    end
    if haskey(basis.extra_data, m)
        return basis.extra_data[m]
    else
        return nothing
    end
end

function normalize(m::DeformationMap{T}) where {T <: SmashProductLieElem}
    nz_index = findfirst(x -> !iszero(x), m)
    if nz_index === nothing
        return m
    end
    lc = leading_coefficient(data(m[CartesianIndex(nz_index[2], nz_index[1])]))
    m = map(e -> 1 // lc * e, m)
end

function filter_independent(R, input)
    deduplicated = unique(Iterators.filter(!iszero, input))::Vector{eltype(input)}
    _, rels = is_linearly_independent_with_relations(R, deduplicated)
    inds = [findlast(!iszero, vec(rels[i, :]))::Int for i in 1:nrows(rels)]
    deleteat!(deduplicated, inds)
    return deduplicated, length(deduplicated)
end
