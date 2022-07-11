mutable struct FreeAlgebra{C <: RingElement} <: Algebra{C}
    base_ring::Ring
    S::Vector{Symbol}
    num_gens::Int

    function FreeAlgebra{C}(R::Ring, symb::Vector{Symbol}) where {C <: RingElement}
        return new{C}(R, symb, length(symb))
    end

    function FreeAlgebra{C}(R::Ring, symb::Vector{String}) where {C <: RingElement}
        return new{C}(R, map(Symbol, symb), length(symb))
    end

end

mutable struct FreeAlgebraElem{C <: RingElement} <: AlgebraElem{C}
    coeffs::Vector{C}
    monoms::Vector{Vector{Int}}
    length::Int
    parent::FreeAlgebra{C}

    function FreeAlgebraElem{C}(A::FreeAlgebra{C}) where {C <: RingElement}
        return new{C}(Array{C}(undef, 0), Array{Vector{Int}}(undef, 0), 0, A)
    end

    function FreeAlgebraElem{C}(A::FreeAlgebra{C}, c::Vector{C}, m::Vector{Vector{Int}}) where {C <: RingElement}
        length(c) == length(m) || throw(DimensionMismatch("c and m are requiered to have the same length."))
        zeroinds = findall(iszero, c)
        deleteat!(c, zeroinds)
        deleteat!(m, zeroinds)
        return new{C}(c, m, length(c), A)
    end

    function FreeAlgebraElem{C}(A::FreeAlgebra{C}, a::C) where {C <: RingElement}
        return new{C}([a], [Int[]], 1, A)
    end

end

function free_algebra(R::Ring, S::Vector{Symbol})
    alg = FreeAlgebra{elem_type(R)}(R, S)
    return alg, gens(alg)
end

function free_algebra(R::Ring, S::Vector{String})
    alg = FreeAlgebra{elem_type(R)}(R, S)
    return alg, gens(alg)
end

parent_type(::Type{FreeAlgebraElem{C}}) where {C <: RingElement} = FreeAlgebra{C}

elem_type(::Type{FreeAlgebra{C}}) where {C <: RingElement} = FreeAlgebraElem{C}


function is_gen(a::FreeAlgebraElem{C}) where {C <: RingElement}
    return length(a) == 1 && isone(a.coeffs[1]) && length(a.monoms[1]) == 1
end

function is_monomial(a::FreeAlgebraElem{C}) where {C <: RingElement}
    return length(a) == 1 && isone(a.coeffs[1])
end

function iszero(a::FreeAlgebraElem)
    return a.length == 0
end

function isone(a::FreeAlgebraElem)
    return a.length == 1 && isempty(a.monoms[1]) && isone(a.coeffs[1])
end

function Base.deepcopy_internal(a::FreeAlgebraElem{C}, dict::IdDict) where {C <: RingElement}
    Rm = deepcopy_internal(a.monoms, dict)
    Rc = Array{C}(undef, a.length)
    for i in 1:a.length
        Rc[i] = deepcopy(a.coeffs[i])
    end
    return parent(a)(Rc, Rm)
end

function showname(::Type{FreeAlgebra{C}}) where {C <: RingElement}
    return "Free Algebra"
end


function Base.:(==)(A1::FreeAlgebra{C}, A2::FreeAlgebra{C}) where {C <: RingElement}
    return (A1.base_ring, A1.S, A1.num_gens) == (A2.base_ring, A2.S, A2.num_gens)
end


function change_base_ring(R::Ring, A::FreeAlgebra{C}) where {C <: RingElement}
    return free_algebra(R, A.S)[1]
end
