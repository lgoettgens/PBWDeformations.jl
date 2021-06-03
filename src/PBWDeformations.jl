module PBWDeformations

abstract type SmashProduct end

include("SmashProductsLie.jl")

"""
    add(x:Integer, y:Integer) :: Integer

adds two integer numbers
```jldoctest
julia> PBWDeformations.add(12,30)
42
```
"""
add(x::Integer, y::Integer) = x+y

end
