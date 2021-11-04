module PBWDeformations

using Combinatorics
using Oscar

import AbstractAlgebra: NCRing, NCRingElem, Ring, RingElement, base_ring, check_parent, coeff, elem_type, gen, gens, isgen, ismonomial, monomial, ngens, parent_type, symbols

import Base: Array, deepcopy, deepcopy_internal, hash, isone, iszero, length, one, parent, show, xor, zero, +, -, *, ^, ==

GAP = Oscar.GAP.Globals
toGAP = Oscar.GAP.julia_to_gap
fromGAP = Oscar.GAP.gap_to_julia

include("Util.jl")
include("Algebra.jl")
include("FreeAlgebra.jl")
include("QuadraticQuoAlgebra.jl")


function sanitize_lie_input(dynkin::Char, n::Int64) :: Nothing
    @assert dynkin in ['A', 'B', 'C', 'D']
    if dynkin == 'A'
        @assert n >= 1
    elseif dynkin == 'B' || dynkin == 'C'
        @assert n >= 2
    elseif dynkin == 'D'
        @assert n >= 4
    end
end


end
