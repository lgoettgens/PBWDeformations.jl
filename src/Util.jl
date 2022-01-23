export flatten, groupBy

"""
    flatten(a::Vector{Vector{T}}) where T
Returns a vector of all elements of elements of `a`.
```jldoctest; setup = :(flatten = PBWDeformations.flatten)
julia> flatten([[1],[],[2,3,4],[5],[]])
5-element Vector{Any}:
 1
 2
 3
 4
 5
```
"""
function flatten(a::Vector{Vector{T}}) where T
    return vcat(a...)
end

"""
    groupBy(a::Vector{T}; eq=(==)) where T
Returns a vector containing the elements of `a` grouped into subvectors of consecutive equal elements.
```jldoctest; setup = :(groupBy = PBWDeformations.groupBy)
julia> groupBy([1,1,2,2,2,2,3,1,4,4])
5-element Vector{Vector{Int64}}:
 [1, 1]
 [2, 2, 2, 2]
 [3]
 [1]
 [4, 4]

julia> groupBy([i for i in -5:5]; eq=((x, y) -> sign(x) == sign(y)))
3-element Vector{Vector{Int64}}:
 [-5, -4, -3, -2, -1]
 [0]
 [1, 2, 3, 4, 5]
```
"""
function groupBy(a::Vector{T}; eq=(==)) where T
    if isempty(a)
        return Vector{T}[]
    end
    v = first(a)
    r = [[v]]
    for x in a[2:end]
        if eq(v, x)
            push!(r[end], x)
        else
            push!(r, [x])
        end
        v = x
    end
    return r
end

function isvalidlie_for_gap(dynkin::Char, n::Int64)
    if dynkin == 'A'
        return n >= 1
    elseif dynkin == 'B'
        return n >= 2
    elseif dynkin == 'C'
        return n >= 2
    elseif dynkin == 'D'
        return n >= 4
    elseif dynkin == 'E'
        return 6 <= n <= 8
    elseif dynkin == 'F'
        return n == 4
    elseif dynkin == 'G'
        return n == 2
    else
        return false
    end
end
