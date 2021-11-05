
mutable struct QuadraticQuoAlgebra{C <: RingElement} <: Algebra{C}
    base_ring :: Ring
    S :: Vector{Symbol}
    num_gens :: Int
    free_alg :: FreeAlgebra{C}
    rels #:: Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}

    function QuadraticQuoAlgebra{C}(free_alg::FreeAlgebra{C}, rels::Dict{Tuple{Int,Int},FreeAlgebraElem{C}}) where C <: RingElement
        this = new{C}(free_alg.base_ring, free_alg.S, free_alg.num_gens, free_alg, Dict{Tuple{Int,Int},FreeAlgebraElem{C}}())
        this.rels = Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}(k => this(a) for (k, a) in rels)
        return this
    end

end

mutable struct QuadraticQuoAlgebraElem{C <: RingElement} <: AlgebraElem{C}
    coeffs :: Vector{C}
    monoms :: Vector{Vector{Int}}
    length :: Int
    parent :: QuadraticQuoAlgebra{C}

    function QuadraticQuoAlgebraElem{C}(A::QuadraticQuoAlgebra) where C <: RingElement
        return new{C}(Array{C}(undef, 0), Array{Vector{Int}}(undef, 0), 0, A)
    end

    function QuadraticQuoAlgebraElem{C}(A::QuadraticQuoAlgebra, c::Vector{C}, m::Vector{Vector{Int}}) where C <: RingElement
        length(c) == length(m) || throw(DimensionMismatch("c and m are requiered to have the same length."))
        zeroinds = findall(iszero, c)
        deleteat!(c, zeroinds)
        deleteat!(m, zeroinds)
        return new{C}(c, m, length(c), A)
    end

    function QuadraticQuoAlgebraElem{C}(A::QuadraticQuoAlgebra, a::C) where C <: RingElement
        return new{C}([a], [Int[]], 1, A)
    end

    function QuadraticQuoAlgebraElem{C}(A::QuadraticQuoAlgebra{C}, b::FreeAlgebraElem{C}) where C <: RingElement
        A.free_alg === parent(b) || throw(ArgumentError("Non-matching algebras"))
        return new{C}(deepcopy(b.coeffs), deepcopy(b.monoms), length(b.coeffs), A)
    end

end

quadratic_quo_algebra(free_alg::FreeAlgebra{C}, rels::Dict{Tuple{Int,Int},FreeAlgebraElem{C}}) where C <: RingElement = QuadraticQuoAlgebra{C}(free_alg, rels)
quo(free_alg::FreeAlgebra{C}, rels::Dict{Tuple{Int,Int},FreeAlgebraElem{C}}) where C <: RingElement = QuadraticQuoAlgebra{C}(free_alg, rels)


parent_type(::Type{QuadraticQuoAlgebraElem{C}}) where C <: RingElement = QuadraticQuoAlgebra{C}

elem_type(::Type{QuadraticQuoAlgebra{C}}) where C <: RingElement = QuadraticQuoAlgebraElem{C}


function isgen(a::QuadraticQuoAlgebraElem{C}) where C <: RingElement
    return length(a) == 1 && isone(a.coeffs[1]) && length(a.monoms[1]) == 1
end

function ismonomial(a::QuadraticQuoAlgebraElem{C}) where C <: RingElement # TODO: fix with normal_form
    return length(a) == 1 && isone(a.coeffs[1])
end

function iszero(a::QuadraticQuoAlgebraElem)
    return a.length == 0 || normal_form(a).length == 0
end

function isone(a::QuadraticQuoAlgebraElem) # TODO: fix with normal_form
    return a.length == 1 && isempty(a.monoms[1]) && isone(a.coeffs[1])
end

function Base.deepcopy_internal(a::QuadraticQuoAlgebraElem{C}, dict::IdDict) where C <: RingElement
    Rm = deepcopy_internal(a.monoms, dict)
    Rc = Array{C}(undef, a.length)
    for i = 1:a.length
       Rc[i] = deepcopy(a.coeffs[i])
    end
    return parent(a)(Rc, Rm)
end

function showname(::Type{QuadraticQuoAlgebra{C}}) where C <: RingElement
    return "Quadratic Quotient Algebra"
end


function Base.:(==)(A1::QuadraticQuoAlgebra{C}, A2::QuadraticQuoAlgebra{C}) where C <: RingElement
    return (A1.base_ring, A1.S, A1.num_gens, A1.free_rels) == (A2.base_ring, A2.S, A2.num_gens, A2.free_rels)
end

function (A::QuadraticQuoAlgebra{C})(b::FreeAlgebraElem{C}) where C <: RingElement
    A.free_alg === parent(b) || throw(ArgumentError("Non-matching algebras"))
    return elem_type(A)(A, b)
end

###############################################################################
#
#   QuadraticQuoAlgebra specific functions
#
###############################################################################

function normal_form(a::QuadraticQuoAlgebraElem{C}) where C <: RingElement
    # TODO
    return a
end
