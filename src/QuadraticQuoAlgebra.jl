mutable struct QuadraticQuoAlgebra{C <: RingElement} <: Algebra{C}
    base_ring :: Ring
    S :: Vector{Symbol}
    num_gens :: Int
    free_alg :: FreeAlgebra{C}
    rels #:: Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}

    function QuadraticQuoAlgebra{C}(free_alg::FreeAlgebra{C}, rels::Dict{Tuple{Int,Int}, FreeAlgebraElem{C}}) where C <: RingElement
        this = new{C}(free_alg.base_ring, free_alg.S, free_alg.num_gens, free_alg, Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}())
        this.rels = Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}(k => this(a) for (k, a) in rels)
        return this
    end

    function QuadraticQuoAlgebra{C}(paren_alg::QuadraticQuoAlgebra{C}, rels) where C <: RingElement #rels::Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}
        this = new{C}(paren_alg.base_ring, paren_alg.S, paren_alg.num_gens, paren_alg.free_alg, Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}())
        this.rels = Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}(k => this(a) for (k, a) in [collect(paren_alg.rels); collect(rels)])
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

    function QuadraticQuoAlgebraElem{C}(A::QuadraticQuoAlgebra{C}, b::QuadraticQuoAlgebraElem{C}) where C <: RingElement
        A.free_alg === parent(b).free_alg || throw(ArgumentError("Non-matching algebras"))
        return new{C}(deepcopy(b.coeffs), deepcopy(b.monoms), length(b.coeffs), A)
    end

    function QuadraticQuoAlgebraElem{C1}(A::QuadraticQuoAlgebra{C1}, b::QuadraticQuoAlgebraElem{C2}) where {C1 <: RingElement, C2 <: RingElement}
        A.S == parent(b).S || throw(ArgumentError("Non-matching algebras"))
        return new{C1}(map(base_ring(A), deepcopy(b.coeffs)), deepcopy(b.monoms), length(b.coeffs), A)
    end

end

function quadratic_quo_algebra(free_alg::FreeAlgebra{C}, rels::Dict{Tuple{Int,Int},FreeAlgebraElem{C}}) where C <: RingElement
    return QuadraticQuoAlgebra{C}(free_alg, rels)
end
function quadratic_quo_algebra(free_alg::QuadraticQuoAlgebra{C}, rels::Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}) where C <: RingElement
    return QuadraticQuoAlgebra{C}(free_alg, rels)
end
function quo(free_alg::FreeAlgebra{C}, rels::Dict{Tuple{Int,Int},FreeAlgebraElem{C}}) where C <: RingElement
    return QuadraticQuoAlgebra{C}(free_alg, rels)
end
function quo(free_alg::QuadraticQuoAlgebra{C}, rels::Dict{Tuple{Int,Int},FreeAlgebraElem{C}}) where C <: RingElement
    return QuadraticQuoAlgebra{C}(free_alg, rels)
end


parent_type(::Type{QuadraticQuoAlgebraElem{C}}) where C <: RingElement = QuadraticQuoAlgebra{C}

elem_type(::Type{QuadraticQuoAlgebra{C}}) where C <: RingElement = QuadraticQuoAlgebraElem{C}


function isgen(a::QuadraticQuoAlgebraElem{C}) where C <: RingElement
    return length(a) == 1 && isone(a.coeffs[1]) && length(a.monoms[1]) == 1
end

# function ismonomial(a::QuadraticQuoAlgebraElem{C}) where C <: RingElement # TODO: fix with normal_form
#     return length(a) == 1 && isone(a.coeffs[1])
# end

Base.:(==)(a::QuadraticQuoAlgebraElem{C}, b::QuadraticQuoAlgebraElem{C}; strict::Bool=false) where C <: RingElement = isequal(a, b; strict)
    
function Base.isequal(a::QuadraticQuoAlgebraElem{C}, b::QuadraticQuoAlgebraElem{C}; strict::Bool=false) where C <: RingElement
    check_parent(a, b, false) || return false
    return iszero(a - b; strict)
end

Base.:(==)(n::Union{Integer, Rational, AbstractFloat}, a::QuadraticQuoAlgebraElem; strict::Bool=false) = isequal(a, n; strict)
Base.isequal(n::Union{Integer, Rational, AbstractFloat}, a::QuadraticQuoAlgebraElem; strict::Bool=false) = isequal(a, n; strict)
Base.:(==)(a::QuadraticQuoAlgebraElem, n::Union{Integer, Rational, AbstractFloat}; strict::Bool=false) = isequal(a, n; strict)

function Base.isequal(a::QuadraticQuoAlgebraElem, n::Union{Integer, Rational, AbstractFloat}; strict::Bool=false)
    return iszero(a - parent(a)(n); strict)
end

Base.:(==)(n::C, a::QuadraticQuoAlgebraElem{C}; strict::Bool=false) where C <: RingElem = isequal(a, n; strict)
Base.:isequal(n::C, a::QuadraticQuoAlgebraElem{C}; strict::Bool=false) where C <: RingElem = isequal(a, n; strict)
Base.:(==)(a::QuadraticQuoAlgebraElem{C}, n::C; strict::Bool=false) where C <: RingElem = isequal(a, n; strict)

function Base.isequal(a::QuadraticQuoAlgebraElem{C}, n::C; strict::Bool=false) where C <: RingElem
    return iszero(a - parent(a)(n); strict)
end

function iszero(a::QuadraticQuoAlgebraElem; strict::Bool=false)
    return strict ? normal_form(a).length == 0 : a.length == 0
end

function isone(a::QuadraticQuoAlgebraElem; strict::Bool=false)
    if strict
        b = normal_form(a)
    else
        b = a
    end
    return b.length == 1 && isempty(b.monoms[1]) && isone(b.coeffs[1])
end

function Base.deepcopy_internal(a::QuadraticQuoAlgebraElem{C}, dict::IdDict) where C <: RingElement
    Rm = deepcopy(a.monoms)
    Rc = Array{C}(undef, a.length)
    for i = 1:a.length
       Rc[i] = deepcopy(a.coeffs[i])
    end
    return parent(a)(Rc, Rm)
end

function showname(::Type{QuadraticQuoAlgebra{C}}) where C <: RingElement
    return "Quadratic Quotient Algebra"
end


function Base.:(==)(A1::QuadraticQuoAlgebra{C}, A2::QuadraticQuoAlgebra{C}) where C <: RingElement # TODO: fix
    return (A1.base_ring, A1.S, A1.num_gens, A1.free_rels) == (A2.base_ring, A2.S, A2.num_gens, A2.free_rels)
end

function (A::QuadraticQuoAlgebra{C})(b::FreeAlgebraElem{C}) where C <: RingElement
    A.free_alg === parent(b) || throw(ArgumentError("Non-matching algebras"))
    return elem_type(A)(A, b)
end

function (A::QuadraticQuoAlgebra{C})(b::QuadraticQuoAlgebraElem{C}) where C <: RingElement
    A.free_alg === parent(b).free_alg || throw(ArgumentError("Non-matching algebras"))
    return elem_type(A)(A, b)
end

function (A::QuadraticQuoAlgebra{C1})(b::QuadraticQuoAlgebraElem{C2}) where {C1 <: RingElement, C2 <: RingElement}
    A.S == parent(b).S || throw(ArgumentError("Non-matching algebras"))
    return elem_type(A)(A, b)
end


function change_base_ring(R::Ring, A::QuadraticQuoAlgebra{C1}) where C1 <: RingElement
    C2 = elem_type(R)
    free_alg = change_base_ring(R, A.free_alg)

    this = quadratic_quo_algebra(free_alg, Dict{Tuple{Int,Int}, FreeAlgebraElem{C2}}())
    this.rels = Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C2}}(k => this(a) for (k, a) in A.rels)

    return this
end

###############################################################################
#
#   QuadraticQuoAlgebra specific functions
#
###############################################################################

function comm(a::QuadraticQuoAlgebraElem{C}, b::QuadraticQuoAlgebraElem{C}; strict=false) where C <: RingElement
    r = a*b - b*a
    return strict ? normal_form(r) : r
end

function normal_form(a::QuadraticQuoAlgebraElem{C}) where C <: RingElement
    todo = deepcopy(a)
    result = zero(parent(todo))
    rels = parent(a).rels
    R = base_ring(a)
    while todo.length > 0
        c = coeff(todo, 1)
        m = monomial(todo, 1)
        exp = m.monoms[1]
        t = c*m
        todo -= t

        changed = false
        for i in 1:length(exp)-1
            if exp[i] > exp[i+1] && haskey(rels, (exp[i], exp[i+1]))
                changed = true
                todo += c * parent(a)([one(R)], [exp[1:i-1]]) * rels[(exp[i], exp[i+1])] * parent(a)([one(R)], [exp[i+2:end]])
                break
            end
        end
        if !changed
            result += t
        end
    end
    return result
end
