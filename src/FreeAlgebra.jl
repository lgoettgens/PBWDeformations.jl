mutable struct FreeAlgebra{C <: RingElement} <: NCRing
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

mutable struct FreeAlgebraElem{C <: RingElement} <: NCRingElem
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


    function FreeAlgebraElem{C1}(
        A::FreeAlgebra{C1},
        b::FreeAlgebraElem{C2},
    ) where {C1 <: RingElement, C2 <: RingElement}
        A.S == parent(b).S || throw(ArgumentError("Non-matching algebras"))
        return new{C1}(map(base_ring(A), deepcopy(b.coeffs)), deepcopy(b.monoms), length(b.coeffs), A)
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

parent(a::FreeAlgebraElem{C}) where {C <: RingElement} = a.parent

base_ring(A::FreeAlgebra{C}) where {C <: RingElement} = A.base_ring::parent_type(C)

base_ring(a::FreeAlgebraElem{C}) where {C <: RingElement} = base_ring(parent(a))

coefficient_ring(A::FreeAlgebra{C}) where {C <: RingElement} = A.base_ring::parent_type(C)

coefficient_ring(a::FreeAlgebraElem{C}) where {C <: RingElement} = coefficient_ring(parent(a))

parent_type(a::FreeAlgebraElem{C}) where {C <: RingElement} = parent_type(typeof(a))

elem_type(A::FreeAlgebra{C}) where {C <: RingElement} = elem_type(typeof(A))

symbols(A::FreeAlgebra) = A.S

ngens(A::FreeAlgebra) = A.num_gens

function gen(A::FreeAlgebra{C}, i::Int) where {C <: RingElement}
    1 <= i <= A.num_gens || throw(ArgumentError("Invalid generator index `i`."))
    return FreeAlgebraElem{C}(A, [one(A.base_ring)], [[i]])
end

function gens(A::FreeAlgebra{C}) where {C <: RingElement}
    return [gen(A, i) for i in 1:A.num_gens]
end

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

# variables actually occuring in a.
function vars(a::FreeAlgebraElem{C}) where {C <: RingElement}
    return [gen(parent(a), i) for i in var_ids(a)]
end

function var_ids(a::FreeAlgebraElem{C}) where {C <: RingElement}
    return sort!(unique!(flatten(a.monoms)))
end

function check_parent(a1::FreeAlgebraElem{C}, a2::FreeAlgebraElem{C}, throw::Bool=true) where {C <: RingElement}
    b = parent(a1) !== parent(a2)
    b & throw && error("Incompatible quadratic algebra elements")
    return !b
end

function coeff(a::FreeAlgebraElem, i::Int)
    return a.coeffs[i]
end

function coeff(a::FreeAlgebraElem{C}, m::Vector{Int}) where {C <: RingElement}
    i = get_index(a, m)
    return isnothing(i) ? zero(C) : a.coeffs[i]
end

function coeff(a::FreeAlgebraElem{C}, m::FreeAlgebraElem{C}) where {C <: RingElement}
    is_monomial(m) || throw(ArgumentError("`m` needs to be a monomial"))
    return coeff(a, m.monoms[1])
end

function get_index(a::FreeAlgebraElem{C}, m::Vector{Int}) where {C <: RingElement}
    return findfirst(x -> x == m, a.monoms)
end

function monomial(a::FreeAlgebraElem, i::Int)
    R = base_ring(a)
    return parent(a)([one(R)], [a.monoms[i]])
end

function term(a::FreeAlgebraElem, i::Int)
    return parent(a)([deepcopy(a.coeffs[i])], [a.monoms[i]])
end

function length(a::FreeAlgebraElem)
    return a.length
end

function zero(A::FreeAlgebra{C}) where {C <: RingElement}
    return A()
end

function zero(a::FreeAlgebraElem)
    return zero(parent(a))
end

function one(A::FreeAlgebra{C}) where {C <: RingElement}
    return A(one(A.base_ring))
end

function one(a::FreeAlgebraElem)
    return one(parent(a))
end

function canonical_unit(a::FreeAlgebraElem{C}) where {C <: RingElement}
    if iszero(a)
        return one(parent(a))
    end
    ind = argmin(a.monoms[1:length(a)])
    return canonical_unit(a.coeffs[ind])
end


function Base.:(==)(A1::FreeAlgebra{C}, A2::FreeAlgebra{C}) where {C <: RingElement}
    return (A1.base_ring, A1.S, A1.num_gens) == (A2.base_ring, A2.S, A2.num_gens)
end

function Base.hash(a::FreeAlgebraElem{C}, h::UInt) where {C <: RingElement}
    b1 = 0x39fec7cd66dff6b4 % UInt
    b2 = 0x1e9a94f334b676fc % UInt
    h1 = xor(h, b1)
    h2 = xor(h, b2)

    p = sortperm(a.monoms)

    for i in 1:length(a)
        h1 = hash(a.coeffs[p[i]], h1)
        h2 = hash(a.monoms[p[i]], h2)
    end

    return xor(h1, h2)
end

Base.:(==)(a::FreeAlgebraElem{C}, b::FreeAlgebraElem{C}) where {C <: RingElement} = isequal(a, b)

function Base.isequal(a::FreeAlgebraElem{C}, b::FreeAlgebraElem{C}) where {C <: RingElement}
    check_parent(a, b, false) || return false
    return iszero(a - b)
end

Base.:(==)(n::Union{Integer, Rational, AbstractFloat}, a::FreeAlgebraElem) = isequal(a, n)
Base.isequal(n::Union{Integer, Rational, AbstractFloat}, a::FreeAlgebraElem) = isequal(a, n)
Base.:(==)(a::FreeAlgebraElem, n::Union{Integer, Rational, AbstractFloat}) = isequal(a, n)

function Base.isequal(a::FreeAlgebraElem, n::Union{Integer, Rational, AbstractFloat})
    return iszero(a - parent(a)(n))
end

Base.:(==)(n::C, a::FreeAlgebraElem{C}) where {C <: RingElem} = isequal(a, n)
Base.:isequal(n::C, a::FreeAlgebraElem{C}) where {C <: RingElem} = isequal(a, n)
Base.:(==)(a::FreeAlgebraElem{C}, n::C) where {C <: RingElem} = isequal(a, n)

function Base.isequal(a::FreeAlgebraElem{C}, n::C) where {C <: RingElem}
    return iszero(a - parent(a)(n))
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

function show(io::IO, A::FreeAlgebra)
    local max_gens = 5 # largest number of generators to print
    n = A.num_gens
    print(io, showname(typeof(A)) * " with ")
    if n > max_gens
        print(io, A.num_gens)
        print(io, " generators ")
    end
    for i in 1:min(n - 1, max_gens - 1)
        print(io, string(A.S[i]), ", ")
    end
    if n > max_gens
        print(io, "..., ")
    end
    print(io, string(A.S[n]))
    print(io, " over ")
    print(IOContext(io, :compact => true), base_ring(A))
end

function show(io::IO, a::FreeAlgebraElem{C}) where {C <: RingElement}
    if iszero(a)
        print(io, "0")
    else
        S = parent(a).S
        for i in 1:length(a)
            if i > 1
                print(io, " + ")
            end
            if isempty(a.monoms[i])
                print(io, a.coeffs[i])
            else
                if !isone(a.coeffs[i])
                    print(io, a.coeffs[i])
                    print(io, "*")
                end
                g = groupBy(a.monoms[i])
                for (j, s) in enumerate(g)
                    if j > 1
                        print(io, "*")
                    end
                    print(io, S[s[1]])
                    if length(s) > 1
                        print(io, "^")
                        print(io, length(s))
                    end
                end
            end
        end
    end
end


function Base.:-(a::FreeAlgebraElem{C}) where {C <: RingElement}
    r = zero(a)
    fit!(r, length(a))
    for i in 1:length(a)
        r.coeffs[i] = -a.coeffs[i]
        r.monoms = deepcopy(a.monoms)
    end
    return r
end


function Base.:+(a::FreeAlgebraElem{C}, b::FreeAlgebraElem{C}) where {C <: RingElement}
    check_parent(a, b)
    r = deepcopy(a)
    for i in 1:length(b)
        j = get_index(a, b.monoms[i])
        if isnothing(j)
            fit!(r, length(r) + 1)
            r.coeffs[end] = deepcopy(b.coeffs[i])
            r.monoms[end] = deepcopy(b.monoms[i])
        else
            r.coeffs[j] += b.coeffs[i]
        end
    end
    zeroinds = findall(iszero, r.coeffs)
    deleteat!(r.coeffs, zeroinds)
    deleteat!(r.monoms, zeroinds)
    r.length -= length(zeroinds)
    return r
end

function Base.:+(a::FreeAlgebraElem{C1}, b::C2) where {C1 <: RingElement, C2 <: RingElement}
    r = deepcopy(a)
    j = get_index(a, Int[])
    if isnothing(j)
        fit!(r, length(r) + 1)
        r.coeffs[end] = deepcopy(b)
        r.monoms[end] = deepcopy(Int[])
    else
        r.coeffs[j] += b
    end
    zeroinds = findall(iszero, r.coeffs)
    deleteat!(r.coeffs, zeroinds)
    deleteat!(r.monoms, zeroinds)
    r.length -= length(zeroinds)
    return r
end

function Base.:+(a::C1, b::FreeAlgebraElem{C2}) where {C1 <: RingElement, C2 <: RingElement}
    return b + a
end


function Base.:-(a::FreeAlgebraElem{C}, b::FreeAlgebraElem{C}) where {C <: RingElement}
    check_parent(a, b)
    r = deepcopy(a)
    for i in 1:length(b)
        j = get_index(a, b.monoms[i])
        if isnothing(j)
            fit!(r, length(r) + 1)
            r.coeffs[end] = -b.coeffs[i]
            r.monoms[end] = deepcopy(b.monoms[i])
        else
            r.coeffs[j] -= b.coeffs[i]
        end
    end
    zeroinds = findall(iszero, r.coeffs)
    deleteat!(r.coeffs, zeroinds)
    deleteat!(r.monoms, zeroinds)
    r.length -= length(zeroinds)
    return r
end


function Base.:*(c::C, a::FreeAlgebraElem{C}) where {C <: RingElem}
    parent(c) === base_ring(a) || throw(ArgumentError("Incompatible coefficient rings"))
    r = zero(a)
    if !iszero(c)
        fit!(r, length(a))
        for i in 1:length(a)
            r.coeffs[i] = c * a.coeffs[i]
            r.monoms = deepcopy(a.monoms)
        end
        zeroinds = findall(iszero, r.coeffs)
        deleteat!(r.coeffs, zeroinds)
        deleteat!(r.monoms, zeroinds)
        r.length -= length(zeroinds)
    end
    return r
end

function Base.:*(c::C, a::FreeAlgebraElem) where {C <: RingElement}
    r = zero(a)
    if !iszero(c)
        fit!(r, length(a))
        for i in 1:length(a)
            r.coeffs[i] = c * a.coeffs[i]
            r.monoms = deepcopy(a.monoms)
        end
        zeroinds = findall(iszero, r.coeffs)
        deleteat!(r.coeffs, zeroinds)
        deleteat!(r.monoms, zeroinds)
        r.length -= length(zeroinds)
    end
    return r
end

function Base.:*(a::FreeAlgebraElem{C}, b::FreeAlgebraElem{C}) where {C <: RingElement}
    check_parent(a, b)
    r = zero(a)
    for i in 1:length(a), j in 1:length(b)
        c = a.coeffs[i] * b.coeffs[j]
        if !iszero(c)
            m = [a.monoms[i]; b.monoms[j]]

            k = get_index(r, m)
            if isnothing(k)
                fit!(r, length(r) + 1)
                r.coeffs[end] = c
                r.monoms[end] = m
            else
                r.coeffs[k] += c
            end
        end
    end
    zeroinds = findall(iszero, r.coeffs)
    deleteat!(r.coeffs, zeroinds)
    deleteat!(r.monoms, zeroinds)
    r.length -= length(zeroinds)
    return r
end


function Base.:^(a::FreeAlgebraElem{C}, n::Int) where {C <: RingElement}
    n >= 0 || throw(DomainError(n, "The exponent needs to be nonnegative."))

    r = one(a)
    while n > 0
        if !iszero(n & 1)
            r = r * a
        end
        n >>= 1
        a = a * a
    end
    zeroinds = findall(iszero, r.coeffs)
    deleteat!(r.coeffs, zeroinds)
    deleteat!(r.monoms, zeroinds)
    r.length -= length(zeroinds)
    return r
end


function comm(a::FreeAlgebraElem{C}, b::FreeAlgebraElem{C}) where {C <: RingElement}
    return a * b - b * a
end

function divexact(a::FreeAlgebraElem{C}, c::C) where {C <: RingElement}
    a = deepcopy(a)
    for i in 1:length(a)
        a.coeffs[i] = divexact(a.coeffs[i], c)
    end
    a
end


function fit!(a::FreeAlgebraElem{C}, n::Int) where {C <: RingElement}
    if length(a) < n
        resize!(a.coeffs, n)
        resize!(a.monoms, n)
        a.length = n
    end
    return nothing
end


function (A::FreeAlgebra{C})(c::Vector{C}, m::Vector{Vector{Int}}) where {C <: RingElement}
    return FreeAlgebraElem{C}(A, c, m)
end

function (A::FreeAlgebra{C})(b::RingElement) where {C <: RingElement}
    return A(base_ring(A)(b))
end

function (A::FreeAlgebra{C})() where {C <: RingElement}
    return FreeAlgebraElem{C}(A)
end

function (A::FreeAlgebra{C})(b::Union{Integer, Rational, AbstractFloat}) where {C <: RingElement}
    return iszero(b) ? A() : FreeAlgebraElem{C}(A, base_ring(A)(b))
end

function (A::FreeAlgebra{C})(b::C) where {C <: RingElement}
    parent(b) === base_ring(A) || throw(ArgumentError("Non-matching base rings"))
    return FreeAlgebraElem{C}(A, b)
end

function (A::FreeAlgebra{C})(b::FreeAlgebraElem{C}) where {C <: RingElement}
    parent(b) === A || throw(ArgumentError("Non-matching algebras"))
    return b
end

function (A::FreeAlgebra{C1})(b::FreeAlgebraElem{C2}) where {C1 <: RingElement, C2 <: RingElement}
    A.S == parent(b).S || throw(ArgumentError("Non-matching algebras"))
    return FreeAlgebraElem{C1}(A, b)
end

function change_base_ring(R::Ring, A::FreeAlgebra{C}) where {C <: RingElement}
    return free_algebra(R, A.S)[1]
end
