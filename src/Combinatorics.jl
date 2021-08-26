function multichoose(n::Integer, k:: Integer) :: Integer
    binomial(n + k - 1, k)
end


# the Multicombinations iterator
struct Multicombinations
    n::Int
    k::Int
end

function Base.iterate(c::Multicombinations, s = [min(1, c.k - i) for i in 1:c.k])
    if c.k < 0 # special case to generate no result for k<0
        return
    end
    if c.k == 0 # special case to generate 1 result for k==0
        isempty(s) && return (s, [1])
        return
    end
    if s[1] == s[c.k] == c.n # end iteration because there is none left
        return
    end
    for i in c.k:-1:1
        if s[i] >= c.n
            continue
        end
        s[i] += 1
        for j in i+1:c.k
            s[j] = s[i]
        end
        break
    end
    return (s, s)
end

Base.length(m::Multicombinations) = multichoose(m.n, m.k)

Base.eltype(::Type{Multicombinations}) = Vector{Int}

"""
    multicombinations(a, k)
Generate all multicombinations of `k` elements from an indexable object `a`. Because the number
of multicombinations can be very large, this function returns an iterator object.
Use `collect(multicombinations(a, k))` to get an array of all multicombinations.
"""
function multicombinations(a, k::Integer)
    reorder(c) = [a[ci] for ci in c]
    (reorder(c) for c in Multicombinations(length(a), k))
end
