module PBWDeformations

using Combinatorics
using Oscar
using SparseArrays

import AbstractAlgebra: NCRing, NCRingElem, Ring, RingElement, base_ring, change_base_ring, check_parent, coeff, elem_type, gen, gens, isgen, ismonomial, monomial, ngens, parent_type, quo, symbols, vars
import Oscar: comm, normal_form

import Base: Array, deepcopy, deepcopy_internal, hash, isone, iszero, length, one, parent, show, xor, zero, +, -, *, ^, ==

GAP = Oscar.GAP

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

include("SmashProductLie.jl")
include("SmashProductDeformLie.jl")


end
