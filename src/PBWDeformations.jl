module PBWDeformations

using Oscar
using SymPy

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

fromSymPy = N

include("Structs.jl")
include("AlgebraElement.jl")
include("QuadraticAlgebra.jl")

# generates
#lie(i::Int64; C::Type=Rational{Int64}) = BasisElement{C}(:lie, i)
#lie(is::Vector{Int64}; C::Type=Rational{Int64}) = Monomial{C}(map(lie, is))
#lie(is::UnitRange{Int64}; C::Type=Rational{Int64}) = lie(collect(is; C))
#lie(i::Int64, j::Int64, ks::Vararg{Int64}; C::Type=Rational{Int64}) = lie([i; j; collect(ks)], C)
#islie(b::BasisElement) = (b[1] === :lie) :: Bool
# and likewise for :grp, :mod, :test

for name in (:lie, :grp, :mod, :test)
    isname = Symbol("is", name)
    @eval begin
        ($name)(i::Int64; C::Type=Rational{Int64}) = BasisElement{C}(Symbol($name), i)
        ($name)(is::Vector{Int64}; C::Type=Rational{Int64}) = Monomial{C}([($name)(i; C) for i in is])
        ($name)(is::UnitRange{Int64}; C::Type=Rational{Int64}) = ($name)(collect(is); C)
        ($name)(i::Int64, j::Int64, ks::Vararg{Int64}; C::Type=Rational{Int64}) = ($name)([i; j; collect(ks)]; C)
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
