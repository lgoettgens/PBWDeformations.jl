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
