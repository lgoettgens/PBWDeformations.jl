abstract type Wrapper end

unpack(x) = x
# Wrappers are supposed to work like haskell's newtype,
# i.e. a "strong" alias which is its own type.
function unpack(w::Wrapper)
    return getfield(w, fieldnames(typeof(w))[1])
end

Base.:(==)(w1::Wrapper, w2::Wrapper) = unpack(w1) == unpack(w2)

Base.getindex(w::Wrapper, ind)  = getindex(unpack(w), ind)
Base.iterate(w::Wrapper)        = iterate(unpack(w))
Base.iterate(w::Wrapper, state) = iterate(unpack(w), state)
Base.keys(w::Wrapper)           = keys(unpack(w))
Base.lastindex(w::Wrapper)      = lastindex(unpack(w))
Base.length(w::Wrapper)         = length(unpack(w))

Base.vcat(X::T...) where {T <: Wrapper} = vcat(map(unpack, X)...)
