module PBWDeformations

using Base: Bool
using Oscar
using SymPy

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

toSymPy = sympify
fromSymPy = N

include("AlgebraElement.jl")
include("QuadraticAlgebra.jl")

lie(i::Int64) = (:lie, i) :: BasisElement
mod(i::Int64) = (:mod, i) :: BasisElement
grp(i::Int64) = (:grp, i) :: BasisElement
isLie(b::BasisElement) = b[1] == :lie
isMod(b::BasisElement) = b[1] == :mod
isGrp(b::BasisElement) = b[1] == :grp

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
