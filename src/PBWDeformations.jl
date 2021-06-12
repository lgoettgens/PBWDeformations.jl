module PBWDeformations

using Oscar

abstract type AlgebraWithCommutators end

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

BasisElement = Tuple{Symbol, Int64}
Coefficient = Integer

lie(i::Int64) = (:lie, i) :: BasisElement
mod(i::Int64) = (:mod, i) :: BasisElement

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
include("SmashProductDeformationLie.jl")

end
