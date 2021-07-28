module PBWDeformations

using Base: _collect
using Oscar
using SymPy

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

fromSymPy = N

include("AlgebraElement.jl")
include("QuadraticAlgebra.jl")

# generates
#lieInt(i::BasisIndex) = (:lie, i) :: BasisElementInternal
#lieInt(is::Vector{BasisIndex}) = map($nameInt, is) :: MonomialInternal
#lieInt(is::UnitRange{BasisIndex}) = lieInt(collect(is)) :: MonomialInternal
#lieInt(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = map($nameInt, [i; j; collect(ks)]) :: MonomialInternal
#lie(i::BasisIndex) = BasisElement(:lie, i) :: BasisElement
#lie(is::Vector{BasisIndex}) = Monomial(map($name, is)) :: Monomial
#lie(is::UnitRange{BasisIndex}) = lie(collect(is)) :: Monomial
#lie(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = Monomial(map($name, [i; j; collect(ks)])) :: Monomial
#islie(b::BasisElement) = (b[1] === :lie) :: Bool
# and likewise for :grp, :mod, :test

for name in (:lie, :grp, :mod, :test)
    nameInt = Symbol(name, "Int")
    isname = Symbol("is", name)
    @eval begin
        ($nameInt)(i::BasisIndex) = (Symbol($name), i) :: BasisElementInternal
        ($nameInt)(is::Vector{BasisIndex}) = map($nameInt, is) :: MonomialInternal
        ($nameInt)(is::UnitRange{BasisIndex}) = ($nameInt)(collect(is)) :: MonomialInternal
        ($nameInt)(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = map($nameInt, [i; j; collect(ks)]) :: MonomialInternal

        ($name)(i::BasisIndex) = BasisElement(Symbol($name), i) :: BasisElement
        ($name)(is::Vector{BasisIndex}) = Monomial(map($name, is)) :: Monomial
        ($name)(is::UnitRange{BasisIndex}) = ($name)(collect(is)) :: Monomial
        ($name)(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = Monomial(map($name, [i; j; collect(ks)])) :: Monomial
        ($isname)(b::BasisElement) = (b[1] === Symbol($name)) :: Bool
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
