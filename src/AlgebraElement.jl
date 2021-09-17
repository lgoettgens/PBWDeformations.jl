abstract type AbstractAlgebraElement <: Wrapper end

ScalarTypes = RingElement
DefaultScalarType = fmpq

struct BasisElement{C <: ScalarTypes} <: AbstractAlgebraElement
    b :: Tuple{Symbol, Int64}

    BasisElement{C}(sym::Symbol, ind::Int64) where C <: ScalarTypes = new{C}((sym, ind))
    BasisElement{C}(b::Tuple{Symbol, Int64}) where C <: ScalarTypes = new{C}(b)
    BasisElement{C}(b::BasisElement{C}) where C <: ScalarTypes = b
end

struct Monomial{C <: ScalarTypes} <: AbstractAlgebraElement
    m :: Vector{BasisElement{C}}

    Monomial{C}() where C <: ScalarTypes = new{C}([])
    Monomial{C}(b::BasisElement{C}) where C <: ScalarTypes = new{C}([b])

    Monomial{C}(m::Vector{BasisElement{C}}) where C <: ScalarTypes = new{C}(m)
    Monomial{C}(m1::BasisElement{C}, ms::Vararg{BasisElement{C}}) where C <: ScalarTypes = new{C}([m1; collect(ms)])
    Monomial{C}(m::Monomial{C}) where C <: ScalarTypes = m

    function Monomial{C}(iterable) :: Monomial{C} where C <: ScalarTypes
        if isempty(iterable)
            return Monomial{C}()
        else
            throw(DomainError("Cannot cast input to Monomial"))
        end
    end
end

struct AlgebraElement{C <: ScalarTypes} <: AbstractAlgebraElement
    a :: Vector{Tuple{C, Monomial{C}}}

    # emtpy list represents empty sum, which is 0
    AlgebraElement{C}() where C <: ScalarTypes = new{C}([])

    AlgebraElement{C}(c::Int64) where C <: ScalarTypes = iszero(c) ? AlgebraElement{C}() : new{C}([(C(c), Monomial{C}())])
    AlgebraElement{C}(c::C) where C <: ScalarTypes = iszero(c) ? AlgebraElement{C}() : new{C}([(c, Monomial{C}())])

    AlgebraElement{C}(b::BasisElement{C}) where C <: ScalarTypes = new{C}([(one(C), Monomial{C}(b))])
    # AlgebraElement{C}(b::BasisElement{C}, one::C) where C <: ScalarTypes = new{C}([(one, Monomial{C}(b))])

    AlgebraElement{C}(m::Monomial{C}) where C <: ScalarTypes = new{C}([(one(C), m)])
    # AlgebraElement{C}(m::Monomial{C}, one::C) where C <: ScalarTypes = new{C}([(one, m)])

    AlgebraElement{C}(a::Vector{Tuple{D, Monomial{C}}}) where {C <: ScalarTypes, D <: C} = new{C}(a)
    AlgebraElement{C}(a::AlgebraElement{C}) where C <: ScalarTypes = a

    function AlgebraElement{C}(iterable) :: AlgebraElement{C} where C <: ScalarTypes
        if isempty(iterable)
            return AlgebraElement{C}()
        else
            throw(DomainError("Cannot cast input to AlgebraElement"))
        end
    end
end

StandardOperand{C <: ScalarTypes} = Union{BasisElement{C}, Monomial{C}, AlgebraElement{C}}
Operand{C <: ScalarTypes} = Union{Int64, <:C, BasisElement{C}, Monomial{C}, AlgebraElement{C}}

function Base.one(::Union{Monomial{C}, Type{Monomial{C}}}) :: Monomial{C} where C <: ScalarTypes
    return Monomial{C}()  
end

function Base.zero(::Union{AlgebraElement{C}, Type{AlgebraElement{C}}}) :: AlgebraElement{C} where C <: ScalarTypes
    return AlgebraElement{C}(0)
end

function Base.one(::Union{AlgebraElement{C}, Type{AlgebraElement{C}}}) :: AlgebraElement{C} where C <: ScalarTypes
    return AlgebraElement{C}(1)
end

function monomials(a::AlgebraElement{C}) :: Vector{Monomial{C}} where C <: ScalarTypes
    return unique([mon for (coeff, mon) in a])
end

function basis_elements(a::AlgebraElement{C}) :: Vector{BasisElement{C}} where C <: ScalarTypes
    return unique(vcat(monomials(a)...))
end


function collect_summands(a::AlgebraElement{C}) :: AlgebraElement{C} where C <: ScalarTypes
    res = Tuple{C, Monomial{C}}[]

    for (coeff, mon) in AlgebraElement{C}(a)
        i = findfirst(x -> x[2] == mon, res)

        if i === nothing
            push!(res, (coeff, mon))
        else
            res[i] = (res[i][1] + coeff, res[i][2])
        end
    end

    return AlgebraElement{C}(filter!(x -> !iszero(x[1]), res))
end

function collect_summands(as::Vector{AlgebraElement{C}}) :: AlgebraElement{C} where C <: ScalarTypes
    return collect_summands(AlgebraElement{C}(vcat(as...)))
end


function issamesum(a1::Operand{C}, a2::Operand{C}) :: Bool where C <: ScalarTypes
    return issetequal(collect_summands(AlgebraElement{C}(a1)),
                      collect_summands(AlgebraElement{C}(a2)))
end

≐ = issamesum # type symbol via \doteq, autocomplete with tab

function Base.iszero(a::AlgebraElement{C}) :: Bool where C <: ScalarTypes
    return a ≐ 0
end

function Base.isone(a::AlgebraElement{C}) :: Bool where C <: ScalarTypes
    return a ≐ 1
end


function Base.:(+)(x::Operand{C}, y::StandardOperand{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return collect_summands([AlgebraElement{C}(x), AlgebraElement{C}(y)])
end

function Base.:(+)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C <: ScalarTypes
    return c + a
end


function Base.:(*)(x::BasisElement{C}, y::BasisElement{C}) :: Monomial{C} where C <: ScalarTypes
    return Monomial{C}([x, y])
end

function Base.:(*)(x::BasisElement{C}, y::Monomial{C}) :: Monomial{C} where C <: ScalarTypes
    return Monomial{C}([x; unpack(y)])
end

function Base.:(*)(x::Monomial{C}, y::BasisElement{C}) :: Monomial{C} where C <: ScalarTypes
    return Monomial{C}([unpack(x); y])
end

function Base.:(*)(x::Monomial{C}, y::Monomial{C}) :: Monomial{C} where C <: ScalarTypes
    return Monomial{C}([unpack(x); unpack(y)])
end

function Base.:(*)(c::C, b::BasisElement{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return AlgebraElement{C}([(c, Monomial{C}(b))])
end

function Base.:(*)(c::C, m::Monomial{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return AlgebraElement{C}([(c, m)])
end

function Base.:(*)(c::Union{Int64, C}, a::AlgebraElement{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return AlgebraElement{C}([(c*coeff, mon) for (coeff, mon) in AlgebraElement{C}(a)])
end

function Base.:(*)(c::Int64, b::BasisElement{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return C(c) * b
end

function Base.:(*)(c::Int64, m::Monomial{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return C(c) * m
end

function Base.:(*)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C <: ScalarTypes
    return c * a
end

function Base.:(*)(a1::AlgebraElement{C}, a2::AlgebraElement{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return collect_summands(AlgebraElement{C}([(c1*c2, m1*m2) for (c1, m1) in a1 for (c2, m2) in a2]))
end

function Base.:(*)(a::AlgebraElement{C}, m::Union{BasisElement{C}, Monomial{C}}) :: AlgebraElement{C} where C <: ScalarTypes
    return collect_summands(AlgebraElement{C}([(c, m2*m) for (c, m2) in a]))
end

function Base.:(*)(m::Union{BasisElement{C}, Monomial{C}}, a::AlgebraElement{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return collect_summands(AlgebraElement{C}([(c, m*m2) for (c, m2) in a]))
end


function Base.:(^)(m::Union{BasisElement{C}, Monomial{C}}, n::Int64) :: Monomial{C} where C <: ScalarTypes
    @assert n >= 0
    return prod([m for _ in 1:n]; init=one(Monomial{C}))
end

function Base.:(^)(a::AlgebraElement{C}, n::Int64) :: AlgebraElement{C} where C <: ScalarTypes
    @assert n >= 0
    return prod([a for _ in 1:n])
end


function Base.:(-)(a::StandardOperand{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return (-1)*a
end

function Base.:(-)(x::Operand{C}, y::StandardOperand{C}) :: AlgebraElement{C} where C <: ScalarTypes
    return x + (-y)
end

function Base.:(-)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C <: ScalarTypes
    return a + (-c)
end


comm(x, y) = x*y - y*x


function change_c(C2::Type{<:ScalarTypes}, b::BasisElement{C1}) :: BasisElement{C2} where C1 <: ScalarTypes
    return BasisElement{C2}(unpack(b))
end

function change_c(C2::Type{<:ScalarTypes}, m::Monomial{C1}) :: Monomial{C2} where C1 <: ScalarTypes
    return Monomial{C2}(map(b -> change_c(C2, b), unpack(m)))
end


function groupby(pred, v)
    if isempty(v)
        return typeof(v)()
    end

    tmp = [0; findall([!pred(v[i], v[i+1]) for i in 1:length(v)-1]); length(v)]

    return [v[tmp[i]+1:tmp[i+1]] for i in 1:length(tmp)-1]
end

function prettyprint(b::BasisElement{C}) :: String where C <: ScalarTypes
    return string(b[1], '(', b[2], ')')
end

function prettyprint(m::Monomial{C}) :: String where C <: ScalarTypes
    if isone(m)
        return "[]"
    else
        g = groupby((x,y) -> x[1] === y[1], m)
        innerFormat = map(l -> string(l[1][1], '(', join([string(t[2]) for t in l], ", "), ')'), g)
        return join(innerFormat, "")
    end
end

function prettyprint(t::Tuple{C, Monomial{C}}) :: String where C <: ScalarTypes
    if isone(t[2])
        # TODO: isinteger is not supported for MPolyElem
        # return string(isinteger(t[1]) ? Int(t[1]) : t[1])
        return string(t[1])
    elseif isone(t[1])
        return prettyprint(t[2])
    elseif isone(-t[1])
        return string('-', prettyprint(t[2]))
    elseif C <: Int64 || C <: Rational{Int64}
        return string(t[1], '⋅', prettyprint(t[2]))
    else
        return string('(', t[1], ')', '⋅', prettyprint(t[2]))
    end
end

function prettyprint(a::AlgebraElement{C}) :: String where C <: ScalarTypes
    if iszero(a)
        return "0"
    else
        return replace(join(map(prettyprint, a), " + "), "+ -" => "- ")
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::StandardOperand{C}) :: Nothing where C <: ScalarTypes
    print(io, prettyprint(x))
end
