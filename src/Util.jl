export flatten, groupBy, isvaliddynkin

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

"""
    groupBy(a::Vector{T}; eq=(==)) where {T}

Returns a vector containing the elements of `a` grouped into subvectors of consecutive equal elements.

# Examples
```jldoctest
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
function groupBy(a::Vector{T}; eq=(==)) where {T}
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


"""
    isvaliddynkin(dynkin::Char, n::Int)

Returns true, if there given parameters uniquely define a dynkin diagram,
i.e. are of one of the forms
  * ``A_n`` for ``n \\geq 1``,
  * ``B_n`` for ``n \\geq 2``,
  * ``C_n`` for ``n \\geq 2``,
  * ``D_n`` for ``n \\geq 4``,
  * ``E_5``, ``E_6``, ``E_7``,
  * ``F_4``,
  * ``G_2``.

# Examples
```jldoctest
julia> isvaliddynkin('A', 2)
true

julia> isvaliddynkin('F', 4)
true

julia> isvaliddynkin('D', 3)
false
```
"""
function isvaliddynkin(dynkin::Char, n::Int)
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
