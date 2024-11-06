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
