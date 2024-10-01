Base.eltype(::Type{<:DeformBasis{C}}) where {C <: RingElem} = DeformationMap{C}

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
    lc = leading_coefficient(m[CartesianIndex(nz_index[2], nz_index[1])])
    m = map(e -> C(1 // lc) * e, m)
end
