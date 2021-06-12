module PBWDeformations

using Oscar

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

BasisElement = Tuple{Symbol, Int64}
Coefficient = Integer
LinearCombination = Vector{Tuple{Coefficient, BasisElement}}

struct AlgebraWithCommutators{T}
    basis :: Vector{BasisElement}

    """
    Contains the simplified commutator of two elements as a linear combination.
    An empty list thus is the empty sum and means that the two elements commutate.
    An absent entry means that there is only a formal commutator.
    """
    commTable :: Dict{Tuple{BasisElement, BasisElement}, LinearCombination}
    extraData :: T
end

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
include("SmashProductDeformLie.jl")

end
