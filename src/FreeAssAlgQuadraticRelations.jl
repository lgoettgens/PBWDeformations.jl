QuadraticRelations{C} = Dict{Tuple{Int, Int}, FreeAssAlgElem{C}} where {C <: RingElem}

function normal_form(a::FreeAssAlgElem{C}, rels::QuadraticRelations{C}) where {C <: RingElem}
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
