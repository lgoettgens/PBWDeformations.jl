mutable struct FreeAlgebra{C <: RingElement} <: Algebra{C}
    base_ring :: Ring
    S :: Vector{Symbol}
    num_gens :: Int

    function FreeAlgebra{C}(R::Ring, S::Vector{Symbol}) where C <: RingElement
        return new{C}(R, S, length(S))
    end

    function FreeAlgebra{C}(R::Ring, S::Vector{String}) where C <: RingElement
        return new{C}(R, [Symbol(s) for s in S], length(S))
    end

end

mutable struct FreeAlgebraElem{C <: RingElement} <: AlgebraElem{C}
    coeffs :: Vector{C}
    monoms :: Vector{Vector{Int}}
    length :: Int
    parent :: FreeAlgebra{C}

    function FreeAlgebraElem{C}(A::Algebra) where C <: RingElement
        return new{C}(Array{C}(undef, 0), Array{Vector{Int}}(undef, 0), 0, A)
    end

    function FreeAlgebraElem{C}(A::Algebra, c::Vector{C}, m::Vector{Vector{Int}}) where C <: RingElement
        length(c) == length(m) || throw(DimensionMismatch("c and m are requiered to have the same length."))
        return new{C}(c, m, length(c), A)
    end

    function FreeAlgebraElem{C}(A::Algebra, a::C) where C <: RingElement
        return new{C}([a], [Int[]], 1, A)
    end

end


parent_type(::Type{FreeAlgebraElem{C}}) where C <: RingElement = FreeAlgebra{C}

elem_type(::Type{FreeAlgebra{C}}) where C <: RingElement = FreeAlgebraElem{C}


function isgen(a::FreeAlgebraElem{C}) where C <: RingElement
    return length(a) == 1 && isone(a.coeffs[1]) && length(a.monoms[1]) == 1
end

function ismonomial(a::FreeAlgebraElem{C}) where C <: RingElement
    return length(a) == 1 && isone(a.coeffs[1])
end

function iszero(a::FreeAlgebraElem)
    return a.length == 0
end

function isone(a::FreeAlgebraElem)
    return a.length == 1 && isempty(a.monoms[1]) && isone(a.coeffs[1])
end

function Base.deepcopy_internal(a::FreeAlgebraElem{C}, dict::IdDict) where C <: RingElement
    Rm = deepcopy_internal(a.monoms, dict)
    Rc = Array{C}(undef, a.length)
    for i = 1:a.length
       Rc[i] = deepcopy(a.coeffs[i])
    end
    return parent(a)(Rc, Rm)
end

function showname(::Type{FreeAlgebra{C}}) where C <: RingElement
    return "Free Algebra"
end


function Base.:(==)(A1::FreeAlgebra{C}, A2::FreeAlgebra{C}) where C <: RingElement
    return (A1.base_ring, A1.S, A1.num_gens) == (A2.base_ring, A2.S, A2.num_gens)
end
