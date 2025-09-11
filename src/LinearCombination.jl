
function Base.push!(lc::LinearCombination{C,T}, pair::Pair{T,C}) where {C<:RingElem, T}
    e, c = pair.first, pair.second
    if haskey(lc.data, e)
        lc.data[e] += c
        if iszero(lc.data[e])
            delete!(lc.data, e)
        end
    else
        lc.data[e] = c
    end
    return lc
end

function Base.iszero(lc::LinearCombination)
    return isempty(lc.data)
end

function Base.show(io::IO, lc::LinearCombination{C,T}) where {C<:RingElem, T}
    if isempty(lc.data)
        print(io, "0")
        return
    end
    join(io, [is_one(c) ? "{$(e)}" : "$(c)*{$(e)}" for (e, c) in lc.data], " + ")
end

function linear_combination(coeffs::Vector{C}, elems::Vector{T}) where {C<:RingElem, T}
    @req length(coeffs) == length(elems) "Length of coeffs and elems must be the same"
    lc = LinearCombination{C,T}()
    for (c, e) in zip(coeffs, elems)
        if !iszero(c)
            push!(lc, (e => c))
        end
    end
    return lc
end

function linear_combination(elems::Vector{Pair{T, C}}) where {C<:RingElem, T}
    lc = LinearCombination{C,T}()
    for (e, c) in elems
        if !iszero(c)
            push!(lc, (e => c))
        end
    end
    return lc
end

function write_in_basis(F::Field, x::T, V::Vector{T}) where {T}
    is_in_span, rel = is_in_span_with_relation(F, x, V)
    @req is_in_span "Element is not in the span"
    return linear_combination(rel, V)
end
