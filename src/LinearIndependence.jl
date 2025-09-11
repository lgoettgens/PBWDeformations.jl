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

"""
    is_linearly_independent(F::Field, V::Vector{T}) where {T} -> Bool

Checks whether the elements of `V` are linearly independent over `F`.
"""
function is_linearly_independent(F::Field, V::Vector{T}) where {T}
    return is_linearly_independent_with_relations(F, V)[1]
end

function is_linearly_independent_with_relations(V::Vector{T}) where {T}
    M = kernel(_linear_independence_coeff_matrix(V); side=:left)
    return nrows(M) == 0, M
end

"""
    is_linearly_independent_with_relations(F::Field, V::Vector{T}) where {T}

Checks whether the elements of `V` are linearly independent over `F`.

This function returns a tuple `(is_independent, relations)` where `is_independent`
is a boolean indicating whether the elements are linearly independent and `relations`
is a matrix whose rows are the coefficients of the linear relations.
"""
function is_linearly_independent_with_relations(F::Field, V::Vector{T}) where {T}
    M = kernel(_linear_independence_coeff_matrix(F, V); side=:left)
    return nrows(M) == 0, M
end

"""
    is_in_span(F::Field, x::T, V::Vector{T}) where {T} -> Bool

Checks whether `x` is in the `F`-span of `V`.
"""
function is_in_span(F::Field, x::T, V::Vector{T}) where {T}
    cm = _linear_independence_coeff_matrix(F, [x; V])
    return rank(cm[2:end, :]) == rank(cm)
end

"""
    is_in_span_with_relation(F::Field, x::T, V::Vector{T}) where {T}

This function returns a tuple `(is_in_span, relation)` where `is_in_span`
is a boolean indicating whether `x` is in the `F`-span of `V` and `relations`
is a vector with (one possibility of) coefficients of `V` that result in `x`.
"""
function is_in_span_with_relation(F::Field, x::T, V::Vector{T}) where {T}
    cm = _linear_independence_coeff_matrix(F, [x; V])
    ker = kernel(cm; side=:left)
    rref!(ker)
    is_in_span = nrows(ker) > 0 && isone(ker[1, 1])
    if is_in_span
        return is_in_span, -ker[1, 2:end]
    else
        return is_in_span, [zero(F) for _ in 1:length(V)]
    end
end

"""
    is_span_subset(F::Field, V::Vector{T}, W::Vector{T}) where {T} -> Bool

Checks whether the `F`-span of `V` is a subset of the `F`-span of `W`.
"""
function is_span_subset(F::Field, V::Vector{T}, W::Vector{T}) where {T}
    cm = _linear_independence_coeff_matrix(F, [V; W])
    return rank(cm[length(V)+1:end, :]) == rank(cm)
end

"""
    is_span_equal(F::Field, V::Vector{T}, W::Vector{T}) where {T} -> Bool

Checks whether the `F`-span of `V` is equal to the `F`-span of `W`.
"""
function is_span_equal(F::Field, V::Vector{T}, W::Vector{T}) where {T}
    cm = _linear_independence_coeff_matrix(F, [V; W])
    return rank(cm[1:length(V), :]) == rank(cm[length(V)+1:end, :]) == rank(cm)
end

function _linear_independence_coeff_matrix(V::Vector{<:FieldElem})
    @req length(V) > 0 "For empty vectors, the field needs to be specified"
    return _linear_independence_coeff_matrix(parent(V[1]), V)
end

function _linear_independence_coeff_matrix(F::Field, V::Vector{Any})
    @req is_empty(V) "Only empty vectors may have eltype Any"
    return matrix(F, 0, 0, V)
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

function _linear_independence_coeff_matrix(F::Field, V::Vector{<:FreeAssociativeAlgebraElem})
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

function _linear_independence_coeff_matrix(F::Field, V::Vector{<:SmashProductLieElem})
    return _linear_independence_coeff_matrix(F, map(e -> data(simplify(e)), V))::dense_matrix_type(F)
end

function _linear_independence_coeff_matrix(F::Field, V::Vector{<:SmashProductLieDeformElem})
    return _linear_independence_coeff_matrix(F, map(e -> data(simplify(e)), V))::dense_matrix_type(F)
end
