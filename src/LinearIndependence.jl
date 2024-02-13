function is_linearly_independent(F::Field, V::Vector{T}) where {T}
    return is_linearly_independent_with_relations(F, V)[1]
end

function is_linearly_independent(V::Vector{T}) where {T}
    return is_linearly_independent_with_relations(V)[1]
end

function is_linearly_independent_with_relations(F::Field, V::Vector{<:FieldElem})
    M = kernel(matrix(F, length(V), 1, V); side=:left)
    return nrows(M) == 0, M
end

function is_linearly_independent_with_relations(F::Field, V::Vector{<:PolyRingElem})
    n = length(V)
    if n == 0
        return false, zero_matrix(F, 0, 0)
    end
    rels = identity_matrix(F, n)
    maxdeg = maximum(degree, V)
    for i in 1:maxdeg
        rels = basis_intersection(rels, is_linearly_independent_with_relations(F, [coeff(v, i) for v in V])[2])
        iszero(rels) && return true, rels
    end
    return false, rels
end

function is_linearly_independent_with_relations(F::Field, V::Vector{<:MatElem})
    n = length(V)
    if n == 0
        return false, zero_matrix(F, 0, 0)
    end
    @req all(size(v) == size(V[1]) for v in V) "Size mismatch"
    rels = identity_matrix(F, n)
    for i in eachindex(V[1])
        rels = basis_intersection(rels, is_linearly_independent_with_relations(F, [v[i] for v in V])[2])
        iszero(rels) && return true, rels
    end
    return false, rels
end

function is_linearly_independent_with_relations(F::Field, V::Vector{<:FreeAssAlgElem})
    coeff_maps = [Dict{Vector{Int}, elem_type(F)}(zip(exponent_words(v), coefficients(v))) for v in V]
    words = reduce(union, keys.(coeff_maps))
    n = length(V)
    rels = identity_matrix(F, n)
    for word in words
        rels = basis_intersection(
            rels,
            is_linearly_independent_with_relations(F, [get(coeff_map, word, zero(F)) for coeff_map in coeff_maps])[2],
        )
        iszero(rels) && return true, rels
    end
    return false, rels
end

function basis_intersection(V::MatElem{T}, W::MatElem{T}) where {T <: FieldElem}
    # using Zassenhaus
    @req base_ring(V) == base_ring(W) "Incompatible base rings"
    @req ncols(V) == ncols(W) "Size mismatch"
    if nrows(V) == 0 || nrows(W) == 0
        return zero(V, 0, ncols(V))
    end
    if V == W
        return V
    end
    block_mat = vcat(hcat(V, V), hcat(W, zero(W)))
    s = rref!(block_mat)
    r = rank(block_mat[:, 1:ncols(V)])
    return block_mat[r+1:s, ncols(V)+1:end]
end
