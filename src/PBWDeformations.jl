module PBWDeformations

using Oscar
using SymPy

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

fromSymPy = N

include("AlgebraElement.jl")
include("QuadraticAlgebra.jl")

function createBasisFunctions(name::Symbol)
    s = Symbol("is", name)
    return quote
        ($name)(i::PBWDeformations.BasisIndex) = (Symbol($name), i) :: PBWDeformations.BasisElement
        ($name)(is::Vector{PBWDeformations.BasisIndex}) = map($name, is) :: PBWDeformations.Monomial{PBWDeformations.BasisElement}
        ($name)(i::PBWDeformations.BasisIndex, j::PBWDeformations.BasisIndex, ks::Vararg{PBWDeformations.BasisIndex}) = 
            map($name, [i; j; collect(ks)]) :: PBWDeformations.Monomial{PBWDeformations.BasisElement}
        ($s)(b::PBWDeformations.BasisElement) = (b[1] === Symbol($name)) :: Bool
    end
end

eval(createBasisFunctions(:lie))
eval(createBasisFunctions(:grp))
eval(createBasisFunctions(:mod))

#lie(i::BasisIndex) = (:lie, i) :: BasisElement
#grp(i::BasisIndex) = (:grp, i) :: BasisElement
#mod(i::BasisIndex) = (:mod, i) :: BasisElement
#lie(is::Vector{BasisIndex}) = map(lie, collect(is)) :: Monomial{BasisElement}
#grp(is::Vector{BasisIndex}) = map(grp, collect(is)) :: Monomial{BasisElement}
#mod(is::Vector{BasisIndex}) = map(mod, collect(is)) :: Monomial{BasisElement}
#lie(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = map(lie, [i; j; collect(ks)]) :: Monomial{BasisElement}
#grp(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = map(grp, [i; j; collect(ks)]) :: Monomial{BasisElement}
#mod(i::BasisIndex, j::BasisIndex, ks::Vararg{BasisIndex}) = map(mod, [i; j; collect(ks)]) :: Monomial{BasisElement}
#isLie(b::BasisElement) = b[1] === :lie
#isGrp(b::BasisElement) = b[1] === :grp
#isMod(b::BasisElement) = b[1] === :mod

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
