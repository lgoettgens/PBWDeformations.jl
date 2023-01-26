QuadraticRelations{C} = Dict{Tuple{Int, Int}, FreeAlgebraElem{C}} where {C <: RingElement}

function normal_form(a::FreeAlgebraElem{C}, rels::QuadraticRelations) where {C <: RingElement}
    todo = deepcopy(a)
    result = zero(parent(todo))
    R = base_ring(a)
    while todo.length > 0
        c = coeff(todo, 1)
        m = monomial(todo, 1)
        exp = m.monoms[1]
        t = c * m
        todo -= t

        changed = false
        for i in 1:length(exp)-1
            if exp[i] > exp[i+1] && haskey(rels, (exp[i], exp[i+1]))
                changed = true
                todo +=
                    c *
                    parent(a)([one(R)], [exp[1:i-1]]) *
                    rels[(exp[i], exp[i+1])] *
                    parent(a)([one(R)], [exp[i+2:end]])
                break
            end
        end
        if !changed
            result += t
        end
    end
    return result
end
