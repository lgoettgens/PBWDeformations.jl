abstract type QuadraticAlgebra{C} <: NCRing end
abstract type QuadraticAlgebraElem{C} <: NCRingElem end

mutable struct GenericQuadraticAlgebra{C <: RingElement} <: QuadraticAlgebra{C}
    base_ring :: Ring
    S :: Vector{Symbol}
    num_gens :: Int

    function GenericQuadraticAlgebra{C}(R::Ring, S::Vector{Symbol}) where C <: RingElement
        return new{C}(R, S, length(S))
    end

    function GenericQuadraticAlgebra{C}(R::Ring, S::Vector{String}) where C <: RingElement
        return new{C}(R, [Symbol(s) for s in S], length(S))
    end
end

mutable struct GenericQuadraticAlgebraElem{C <: RingElement} <: QuadraticAlgebraElem{C}
    coeffs :: Vector{C}
    monoms :: Vector{Vector{Int}}
    length :: Int
    parent :: QuadraticAlgebra{C}

    function GenericQuadraticAlgebraElem{C}(A::QuadraticAlgebra) where C <: RingElement
        return new{C}(Array{C}(undef, 0), Array{Vector{Int}}(undef, 0), 0, A)
    end

    function GenericQuadraticAlgebraElem{C}(A::QuadraticAlgebra, c::Vector{C}, m::Vector{Vector{Int}}) where C <: RingElement
        length(c) == length(m) || throw(DimensionMismatch("c and m are requiered to have the same length."))
        return new{C}(c, m, length(c), A)
    end

    function GenericQuadraticAlgebraElem{C}(A::QuadraticAlgebra, a::C) where C <: RingElement
        return new{C}([a], [Int[]], 1, A)
    end

end


parent(a::QuadraticAlgebraElem{C}) where C <: RingElement = a.parent

parent_type(::Type{QuadraticAlgebraElem{C}}) where C <: RingElement = QuadraticAlgebra{C}

elem_type(::Type{QuadraticAlgebra{C}}) where C <: RingElement = QuadraticAlgebraElem{C}

base_ring(A::QuadraticAlgebra{C}) where C <: RingElement = A.base_ring::parent_type(C)

base_ring(a::QuadraticAlgebraElem{C}) where C <: RingElement = base_ring(parent(a))

coefficient_ring(A::QuadraticAlgebra{C}) where C <: RingElement = A.base_ring::parent_type(C)

coefficient_ring(a::QuadraticAlgebraElem{C}) where C <: RingElement = coefficient_ring(parent(a))

symbols(A::QuadraticAlgebra) = A.S

ngens(A::QuadraticAlgebra) = A.num_gens

function gen(A::QuadraticAlgebra{C}, i::Int) where C <: RingElement 
    1 <= i <= A.num_gens || throw(ArgumentError("Invalid generator index `i`."))
    return GenericQuadraticAlgebraElem{C}(A, [one(C)], [[i]])
end

function gens(A::QuadraticAlgebra{C}) where C <: RingElement
    return [gen(A, i) for i in 1:A.num_gens]
end

# TODO: variables actually occuring in a.
# https://github.com/Nemocas/AbstractAlgebra.jl/blob/76891dd23ccf5b3b8dc0f6c6b00fb7ad64e6acc4/src/generic/MPoly.jl#L83
# function vars(a::QuadraticAlgebraElem{C}) where C <: RingElement
# end

function check_parent(a1::QuadraticAlgebraElem{C}, a2::QuadraticAlgebraElem{C}, throw::Bool = true) where C <: RingElement
    b = parent(a1) != parent(a2)
    b & throw && error("Incompatible quadratic algebra elements")
    return !b
end

# TODO: monomial vector manipulation

function coeff(a::QuadraticAlgebraElem, i::Int)
    return a.coeffs[i]
end

function coeff(a::QuadraticAlgebraElem{C}, m::Vector{Int}) where C <: RingElement
    i = get_index(a, m)
    return isnothing(i) ? zero(C) : a.coeffs[i]
end

function coeff(a::QuadraticAlgebraElem{C}, m::QuadraticAlgebraElem{C}) where C <: RingElement
    ismonomial(m) || throw(ArgumentError("`m` needs to be a monomial"))
    return coeff(a, m.monoms[1])
end

# TODO: setcoeff!

function get_index(a::QuadraticAlgebraElem{C}, m::Vector{Int}) where C <: RingElement
    return findfirst(x -> x == m, a.monoms)
end

function Base.hash(a::QuadraticAlgebraElem{C}, h::UInt) where C <: RingElement
    b1 = 0x39fec7cd66dff6b4%UInt
    b2 = 0x1e9a94f334b676fc%UInt
    h1 = xor(h, b1)
    h2 = xor(h, b2)

    for i in 1:length(a)
        h1 = hash(a.coeffs[i], h1)
        h2 = hash(a.monoms[i], h2)
    end

    return xor(h1, h2)
end

# TODO function isunit(a::QuadraticAlgebraElem{C}) where C <: RingElement

function isgen(a::QuadraticAlgebraElem{C}) where C <: RingElement
    return length(a) == 1 && isone(a.coeffs[1]) && length(a.monoms[1]) == 1
end

function ismonomial(a::QuadraticAlgebraElem{C}) where C <: RingElement
    return length(a) == 1 && isone(a.coeffs[1])
end

function monomial(a::QuadraticAlgebraElem, i::Int)
    R = base_ring(a)
    return parent(a)([one(R)], [a.monoms[i]])
end

# TODO monomial!

function term(a::QuadraticAlgebraElem, i::Int)
    R = base_ring(a)
    return parent(a)([deepcopy(a.coeffs[i])], [a.monoms[i]])
end

function length(a::QuadraticAlgebraElem)
    return a.length
end

function iszero(a::QuadraticAlgebraElem)
    return a.length == 0
end

function isone(a::QuadraticAlgebraElem)
    return a.length == 1 && isempty(a.monoms[1]) && isone(a.coeffs[1])
end

function zero(A::QuadraticAlgebra{C}) where C <: RingElement
    return A()
end

function zero(a::QuadraticAlgebraElem)
    return zero(parent(a))
end

function one(A::QuadraticAlgebra{C}) where C <: RingElement
    return A(one(C))
end
function one(a::QuadraticAlgebraElem)
    return one(parent(a))
end

function Base.deepcopy_internal(a::QuadraticAlgebraElem{C}, dict::IdDict) where C <: RingElement
    Rm = deepcopy_internal(a.monoms, dict)
    Rc = Array{C}(undef, a.length)
    for i = 1:a.length
       Rc[i] = deepcopy(a.coeffs[i])
    end
    return parent(a)(Rc, Rm)
end

function show(io::IO, A::QuadraticAlgebra)
    local max_gens = 5 # largest number of generators to print
    n = A.num_gens
    print(io, "Quadratic Algebra with ")
    if n > max_gens
       print(io, A.num_gens)
       print(io, " generators ")
    end
    for i = 1:min(n - 1, max_gens - 1)
       print(io, string(A.S[i]), ", ")
    end
    if n > max_gens
       print(io, "..., ")
    end
    print(io, string(A.S[n]))
    print(io, " over ")
    print(IOContext(io, :compact => true), base_ring(A))
end

function show(io::IO, a::QuadraticAlgebraElem{C}) where C <: RingElement
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
                for s in g
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

function Base.:(==)(A1::QuadraticAlgebra{C}, A2::QuadraticAlgebra{C}) where C <: RingElement
    return (A1.base_ring, A1.S, A1.num_gens) == (A2.base_ring, A2.S, A2.num_gens)
end

###############################################################################
#
#   Arithmetic functions
#
###############################################################################

function Base.:-(a::QuadraticAlgebraElem{C}) where C <: RingElement
    r = zero(a)
    fit!(r, length(a))
    for i in 1:length(a)
        r.coeffs[i] = -a.coeffs[i]
        r.monoms = deepcopy(a.monoms)
    end
    return r
end

function Base.:+(a::QuadraticAlgebraElem{C}, b::QuadraticAlgebraElem{C}) where C <: RingElement
    check_parent(a, b)
    r = deepcopy(a)
    for i in 1:length(b)
        j = get_index(a, b.monoms[i])
        if isnothing(j)
            fit!(r, length(r)+1)
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


function Base.:-(a::QuadraticAlgebraElem{C}, b::QuadraticAlgebraElem{C}) where C <: RingElement
    check_parent(a, b)
    r = deepcopy(a)
    for i in 1:length(b)
        j = get_index(a, b.monoms[i])
        if isnothing(j)
            fit!(r, length(r)+1)
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

function Base.:*(c::C, a::QuadraticAlgebraElem{C}) where C <: RingElement
    parent(c) == base_ring(a) || throw(ArgumentError("Incompatible coefficient rings"))
    r = zero(a)
    if !iszero(c)
        fit!(r, length(a))
        for i in 1:length(a)
            r.coeffs[i] = c*a.coeffs[i]
            r.monoms = deepcopy(a.monoms)
        end
    end
    return r
end


function Base.:*(a::QuadraticAlgebraElem{C}, b::QuadraticAlgebraElem{C}) where C <: RingElement
    check_parent(a, b)
    r = zero(a)
    for i in 1:length(a), j in 1:length(b)
        c = a.coeffs[i] * b.coeffs[j]
        if !iszero(c)
            m = [a.monoms[i]; b.monoms[j]]

            k = get_index(r, m)
            if isnothing(k)
                fit!(r, length(r)+1)
                r.coeffs[end] = c
                r.monoms[end] = m
            else
                r.coeffs[k] += c
            end
        end
    end
    return r
end


## TODO: continue

###############################################################################
#
#   Unsafe functions
#
###############################################################################

## TODO: continue

function fit!(a::QuadraticAlgebraElem{C}, n::Int) where C <: RingElement
    if length(a) < n
       resize!(a.coeffs, n)
       resize!(a.monoms, n)
       a.length = n
    end
    return nothing
end
 
## TODO: continue

###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (A::QuadraticAlgebra{C})(c::Vector{C}, m::Vector{Vector{Int}}) where C <: RingElement
    return GenericQuadraticAlgebraElem{C}(A, c, m)
end

function (A::QuadraticAlgebra{C})(b::RingElement) where C <: RingElement
    return A(base_ring(A)(b))
end
     
function (A::QuadraticAlgebra{C})() where C <: RingElement
    return GenericQuadraticAlgebraElem{C}(A)
end
     
function (A::QuadraticAlgebra{C})(b::Union{Integer, Rational, AbstractFloat}) where C <: RingElement
    return GenericQuadraticAlgebraElem{C}(A, base_ring(A)(b))
end
     
function (A::QuadraticAlgebra{C})(b::C) where C <: RingElement
    parent(b) != base_ring(A) && throw(ArgumentError("Non-matching base rings"))
    return GenericQuadraticAlgebraElem{C}(A, b)
end
     
function (A::QuadraticAlgebra{C})(b::QuadraticAlgebraElem{C}) where C <: RingElement
    parent(b) != A && error("Non-matching algebras")
    return b
end

# mutable struct SmashProductLie{C <: RingElement} <: QuadraticAlgebra{C}
#     base_ring :: Ring
#     S :: Vector{Symbol}
#     num_gens :: Int

#     """
#     Stores relations of the form ab = c for basis elements a,b.
#     An empty list represents the empty sum, in which case ab = 0.
#     An absent entry means that there is no relation, so we cannot simplify ab.
#     """
#     relTable :: Dict{Tuple{BasisElement{C}, BasisElement{C}}, AlgebraElement{C}}
#     extraData

#     SmashProductLie{C}(basis, relTable, extraData = nothing) where {C <: ScalarTypes} =
#         new{C}(basis, relTable, extraData)
# end



# function Base.:(==)(alg1::QuadraticAlgebra{C}, alg2::QuadraticAlgebra{C}) :: Bool where C <: ScalarTypes
#     (alg1.basis, alg1.relTable, alg1.extraData) ==
#     (alg2.basis, alg2.relTable, alg2.extraData)
# end

# function Base.show(io::IO, alg::QuadraticAlgebra{C}) :: Nothing where C <: ScalarTypes
#     println(io, "Algebra with quadratic relations of dimension ", length(alg.basis))
#     println(io, "Relation table has ", length(alg.relTable), " entries")
#     println(io, "Coefficients of type ", C)
#     println(io, "Extra data:")
#     print(io, alg.extraData)
# end

# function Base.in(b::BasisElement{C}, alg::QuadraticAlgebra{C}) :: Bool where C <: ScalarTypes
#     return b in alg.basis
# end

# function Base.in(m::Monomial{C}, alg::QuadraticAlgebra{C}) :: Bool where C <: ScalarTypes
#     return all(b -> b in alg, m)
# end

# function Base.in(a::AlgebraElement{C}, alg::QuadraticAlgebra{C}) :: Bool where C <: ScalarTypes
#     return all(m in alg for (c, m) in a)
# end


# function normal_form(alg::QuadraticAlgebra{C}, a::AlgebraElement{C}) :: AlgebraElement{C} where C <: ScalarTypes
#     todo = copy(unpack(a))
#     result = AlgebraElement{C}(0)

#     while !isempty(todo)
#         coeff, mon = pop!(todo)

#         changed = false
#         for i in 1:length(mon)-1
#             if haskey(alg.relTable, (mon[i], mon[i+1]))
#                 changed = true

#                 # TODO: something like this: todo += coeff * (mon[1:i-1] * alg.relTable[(mon[i], mon[i+1])] * mon[i+2:end])
#                 todo = unpack(AlgebraElement{C}(todo) + coeff * (Monomial{C}(mon[1:i-1]) * alg.relTable[(mon[i], mon[i+1])] * Monomial{C}(mon[i+2:end])))

#                 break
#             end
#         end

#         if !changed
#             result += coeff * mon
#         end
#     end

#     return result
# end

# function normal_form(alg::QuadraticAlgebra{C}, m::Union{BasisElement{C}, Monomial{C}}) :: AlgebraElement{C} where C <: ScalarTypes
#     return normal_form(alg, AlgebraElement{C}(m))
# end
