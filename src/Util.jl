function groupBy(a::Vector{T}; cmp=(==)) where T
    if isempty(a)
        return Vector{T}[]
    end
    v = first(a)
    r = [[v]]
    for x in a[2:end]
        if cmp(v, x)
            push!(r[end], x)
        else
            push!(r, [x])
        end
        v = x
    end
    return r
end
