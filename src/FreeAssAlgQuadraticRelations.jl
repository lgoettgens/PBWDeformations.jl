QuadraticRelations{C} = Dict{Tuple{Int, Int}, FreeAssAlgElem{C}} where {C <: RingElement}

function normal_form(a::FreeAssAlgElem{C}, rels::QuadraticRelations{C}) where {C <: RingElement}
    todo = deepcopy(a)
    result = zero(parent(todo))
    CR = base_ring(a)
    A = parent(todo)
    while todo.length > 0
        c = leading_coefficient(todo)
        exp = leading_exponent_word(todo)
        t = leading_term(todo)
        todo -= t

        changed = false
        for i in 1:length(exp)-1
            if exp[i] > exp[i+1] && haskey(rels, (exp[i], exp[i+1]))
                changed = true
                todo += A([c], [exp[1:i-1]]) * rels[(exp[i], exp[i+1])] * A([one(CR)], [exp[i+2:end]])
                break
            end
        end
        if !changed
            result += t
        end
    end
    return result
end

function comm(a::FreeAssAlgElem{C}, b::FreeAssAlgElem{C}) where {C <: RingElement}
    return a * b - b * a
end


# can be removed once available in AbstractAlgebra upstream
function canonical_unit(a::FreeAssAlgElem{T}) where {T <: RingElement}
    return canonical_unit(leading_coefficient(a))
end

function _change_freeassalg_ring(R, Rx, cached)
    P, _ = AbstractAlgebra.FreeAssociativeAlgebra(R, symbols(Rx); cached=cached)
    return P
end

function change_base_ring(
    R::Ring,
    p::FreeAssAlgElem{T};
    cached=true,
    parent::AbstractAlgebra.FreeAssAlgebra=_change_freeassalg_ring(R, parent(p), cached),
) where {T <: RingElement}
    base_ring(parent) != R && error("Base rings do not match.")
    return _map(R, p, parent)
end

function _map(g, p::FreeAssAlgElem{T}, Rx) where {T <: RingElement}
    cvzip = zip(coefficients(p), exponent_words(p))
    M = MPolyBuildCtx(Rx)
    for (c, v) in cvzip
        push_term!(M, g(c), v)
    end

    return finish(M)
end
