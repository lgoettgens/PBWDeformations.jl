module PBWDeformations

using Oscar
using SymPy

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

fromSymPy = N

include("AlgebraElement.jl")
include("QuadraticAlgebra.jl")

# generates
# lie(i::BasisIndex) = (:lie, i) :: BasisElement
# lie(is::Vector{BasisIndex}) = map(lie, collect(is)) :: Monomial{BasisElement}
# lie(is::UnitRange{BasisIndex}) = lie(collect(is)) :: Monomial{BasisElement}
# lie(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = map(lie, [i; j; collect(ks)]) :: Monomial{BasisElement}
# islie(b::BasisElement) = b[1] === :lie
# and likewise for :grp, :mod, :test

for name in (:lie, :grp, :mod, :test)
    s = Symbol("is", name)
    @eval begin
        ($name)(i::BasisIndex) = (Symbol($name), i) :: BasisElement
        ($name)(is::Vector{BasisIndex}) = map($name, is) :: Monomial{BasisElement}
        ($name)(is::UnitRange{BasisIndex}) = ($name)(collect(is)) :: Monomial{BasisElement}
        ($name)(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = map($name, [i; j; collect(ks)]) :: Monomial{BasisElement}
        ($s)(b::BasisElement) = (b[1] === Symbol($name)) :: Bool
    end
end

function sanitizeLieInput(dynkin::Char, n::Int64) :: Nothing
    @assert dynkin in ['A', 'B', 'C', 'D']
    if dynkin == 'A'
        @assert n >= 1
    elseif dynkin == 'B' || dynkin == 'C'
        @assert n >= 2
    elseif dynkin == 'D'
        @assert n >= 4
    end
end

include("SmashProductLie.jl")
include("SmashProductDeformLie.jl")
include("GroupAlgebra.jl")

end
