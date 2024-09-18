const _linear_independence_rref_cutoff = 0

function column_rref!(mat::MatElem{T}) where {T <: FieldElem}
    rk = rref!(Oscar.AbstractAlgebra.Solve.lazy_transpose(mat))
    return view(mat, :, 1:rk)
end

function column_rref!(mat::Nemo._MatTypes) # change to Nemo._FieldMatTypes
    trmat = transpose!(mat)
    rk = rref!(trmat)
    return view(transpose!(trmat), :, 1:rk)
end

function is_linearly_independent(V::Vector{T}) where {T}
    return is_linearly_independent_with_relations(V)[1]
end

function is_linearly_independent(F::Field, V::Vector{T}) where {T}
    return is_linearly_independent_with_relations(F, V)[1]
end

function is_linearly_independent_with_relations(V::Vector{T}) where {T}
    M = kernel(_linear_independence_coeff_matrix(V); side=:left)
    return nrows(M) == 0, M
end

function is_linearly_independent_with_relations(F::Field, V::Vector{T}) where {T}
    M = kernel(_linear_independence_coeff_matrix(F, V); side=:left)
    return nrows(M) == 0, M
end

function _linear_independence_coeff_matrix(V::Vector{<:FieldElem})
    @req length(V) > 0 "For empty vectors, the field needs to be specified"
    return _linear_independence_coeff_matrix(parent(V[1]), V)
end

function _linear_independence_coeff_matrix(F::Field, V::Vector{<:FieldElem})
    return matrix(F, length(V), 1, V)
end

function _linear_independence_coeff_matrix(F::Field, V::Vector{<:PolyRingElem})
    @req isempty(V) || all(v -> parent(v) === parent(V[1]), V) "Incompatible polynomial rings"
    return reduce(
        hcat,
        begin
            mat = _linear_independence_coeff_matrix(F, [coeff(v, i) for v in V])
            if ncols(mat) > _linear_independence_rref_cutoff * nrows(mat)
                column_rref!(mat)
            else
                mat
            end
        end for i in 0:maximum(degree, V; init=-1);
        init=zero_matrix(F, length(V), 0),
    )::dense_matrix_type(F)
end

function _linear_independence_coeff_matrix(F::Field, V::Vector{<:MatElem})
    n = length(V)
    if n == 0
        return zero_matrix(F, n, 0)
    end
    @req all(v -> parent(v) === parent(V[1]), V) "Incompatible matrix spaces"
    return reduce(
        hcat,
        begin
            mat = _linear_independence_coeff_matrix(F, [v[i] for v in V])
            if ncols(mat) > _linear_independence_rref_cutoff * nrows(mat)
                column_rref!(mat)
            else
                mat
            end
        end for i in eachindex(V[1]);
        init=zero_matrix(F, n, 0),
    )::dense_matrix_type(F)
end

function _linear_independence_coeff_matrix(F::Field, V::Vector{<:FreeAssAlgElem})
    n = length(V)
    if n == 0
        return zero_matrix(F, n, 0)
    end
    R = parent(V[1])
    C = base_ring(R)
    @req all(v -> parent(v) === R, V) "Incompatible algebras"
    coeff_maps = [Dict{Vector{Int}, elem_type(F)}(zip(exponent_words(v), coefficients(v))) for v in V]
    support_words = reduce(union, keys.(coeff_maps))
    return reduce(
        hcat,
        begin
            mat = _linear_independence_coeff_matrix(F, [get(coeff_map, word, zero(C)) for coeff_map in coeff_maps])
            if ncols(mat) > _linear_independence_rref_cutoff * nrows(mat)
                column_rref!(mat)
            else
                mat
            end
        end for word in support_words;
        init=zero_matrix(F, n, 0),
    )::dense_matrix_type(F)
end

