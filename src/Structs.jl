abstract type Wrapper end

# Wrappers are supposed to work like haskell's newtype,
# i.e. a "strong" alias which is its own type.
function unpack(w::Wrapper)
    if iszero(fieldcount(typeof(w)))
        return nothing
    else
        return getfield(w, fieldnames(typeof(w))[1])
    end
end

function Base.:(==)(w1::Wrapper, w2::Wrapper) :: Bool
    T = typeof(w1)
    if typeof(w2) !== T
        return false
    end

    for name in fieldnames(T)
        if getfield(w1, name) != getfield(w2, name)
            return false
        end
    end

    return true
end

function Base.getindex(w::Wrapper, ind)
    return getindex(unpack(w), ind)
end

function Base.iterate(w::Wrapper)
    return iterate(unpack(w))
end

function Base.iterate(w::Wrapper, state)
    return iterate(unpack(w), state)
end

function Base.keys(w::Wrapper)
    return keys(unpack(w))
end

function Base.lastindex(w::Wrapper)
    return lastindex(unpack(w))
end

function Base.length(w::Wrapper)
    return length(unpack(w))
end
