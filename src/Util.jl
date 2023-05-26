"""
    flatten(a::Vector{Vector{T}}) where {T}

Returns a vector of all elements of elements of `a`.

# Example
```jldoctest
julia> flatten([[1],[],[2,3,4],[5],[]])
5-element Vector{Any}:
 1
 2
 3
 4
 5
```
"""
function flatten(a::Vector{Vector{T}}) where {T}
    return vcat(a...)
end

# remove once https://github.com/thofma/Hecke.jl/pull/1085 is available
macro vprintln(s, msg)
    quote
        @vprint $s ($(esc(msg)) * '\n')
    end
end

macro vprintln(s, l::Int, msg)
    quote
        @vprint $s $l ($(esc(msg)) * '\n')
    end
end
