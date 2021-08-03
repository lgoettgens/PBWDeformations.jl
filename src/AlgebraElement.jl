abstract type AbstractAlgebraElement <: Wrapper end

struct BasisElement{C} <: AbstractAlgebraElement
    b :: Tuple{Symbol, Int64}

    BasisElement{C}(sym::Symbol, ind::Int64) where C = new{C}((sym, ind))
    BasisElement{C}(b::Tuple{Symbol, Int64}) where C = new{C}(b)
    BasisElement{C}(b::BasisElement{C}) where C = b
end

struct Monomial{C} <: AbstractAlgebraElement
    m :: Vector{BasisElement{C}}

    Monomial{C}() where C = new{C}([])
    Monomial{C}(b::BasisElement{C}) where C = new{C}([b])

    Monomial{C}(m::Vector{BasisElement{C}}) where C = new{C}(m)
    Monomial{C}(m1::BasisElement{C}, ms::Vararg{BasisElement{C}}) where C = new{C}([m1; collect(ms)])
    Monomial{C}(m::Monomial{C}) where C = m

    function Monomial{C}(a::Vector{<:Any}) :: Monomial{C} where C
        if isempty(a)
            return Monomial{C}()
        else
            throw(DomainError("Cannot cast input to Monomial"))
        end
    end
end

struct AlgebraElement{C} <: AbstractAlgebraElement
    a :: Vector{Tuple{C, Monomial{C}}}

    # emtpy list represents empty sum, which is 0
    AlgebraElement{C}() where C = new{C}([])

    AlgebraElement{C}(c::Int64) where C = iszero(c) ? AlgebraElement{C}() : new{C}([(C(c), Monomial{C}())])
    AlgebraElement{C}(c::C) where C = iszero(c) ? AlgebraElement{C}() : new{C}([(c, Monomial{C}())])

    AlgebraElement{C}(b::BasisElement{C}) where C = new{C}([(C(1), Monomial{C}(b))])
    AlgebraElement{C}(b::BasisElement{C}, one::C) where C = new{C}([(one, Monomial{C}(b))])

    AlgebraElement{C}(m::Monomial{C}) where C = new{C}([(C(1), m)])
    AlgebraElement{C}(m::Monomial{C}, one::C) where C = new{C}([(one, m)])

    AlgebraElement{C}(a::Vector{Tuple{D, Monomial{C}}}) where {C, D <: C} = new{C}(a)
    AlgebraElement{C}(a::AlgebraElement{C}) where C = a

    function AlgebraElement{C}(a::Vector{<:Any}) :: AlgebraElement{C} where C
        if isempty(a)
            return AlgebraElement{C}()
        else
            throw(DomainError("Cannot cast input to AlgebraElement"))
        end
    end
end

StandardOperand{C} = Union{BasisElement{C}, Monomial{C}, AlgebraElement{C}}
Operand{C} = Union{Int64, <:C, BasisElement{C}, Monomial{C}, AlgebraElement{C}}

function Base.one(::Union{Monomial{C}, Type{Monomial{C}}}) :: Monomial{C} where C
    return Monomial{C}()  
end

function Base.zero(::Union{AlgebraElement{C}, Type{AlgebraElement{C}}}) :: AlgebraElement{C} where C
    return AlgebraElement{C}(0)
end

function Base.one(::Union{AlgebraElement{C}, Type{AlgebraElement{C}}}) :: AlgebraElement{C} where C
    return AlgebraElement{C}(1)
end

function monomials(a::AlgebraElement{C}) :: Vector{Monomial{C}} where C
    return unique([mon for (coeff, mon) in a])
end

function basisElements(a::AlgebraElement{C}) :: Vector{BasisElement{C}} where C
    return unique(vcat(monomials(a)...))
end


function collectSummands(a::AlgebraElement{C}) :: AlgebraElement{C} where C
    res = Tuple{C, Monomial{C}}[]

    for (coeff, mon) in AlgebraElement{C}(a)
        i = findfirst(x -> x[2] == mon, res)

        if i === nothing
            push!(res, (coeff, mon))
        else
            res[i] = (res[i][1] + coeff, res[i][2])
        end
    end

    return AlgebraElement{C}(filter(x -> !iszero(x[1]), res))
end

function collectSummands(as::Vector{AlgebraElement{C}}) :: AlgebraElement{C} where C
    return collectSummands(AlgebraElement{C}(vcat(as...)))
end


function sameSum(a1::Operand{C}, a2::Operand{C}) :: Bool where C
    return issetequal(collectSummands(AlgebraElement{C}(a1)),
                      collectSummands(AlgebraElement{C}(a2)))
end

≐ = sameSum # type symbol via \doteq, autocomplete with tab

function Base.iszero(a::AlgebraElement{C}) :: Bool where C
    return a ≐ 0
end

function Base.isone(a::AlgebraElement{C}) :: Bool where C
    return a ≐ 1
end


function Base.:(+)(x::Operand{C}, y::StandardOperand{C}) :: AlgebraElement{C} where C
    return collectSummands([AlgebraElement{C}(x), AlgebraElement{C}(y)])
end

function Base.:(+)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return c + a
end


function Base.:(*)(x::BasisElement{C}, y::BasisElement{C}) :: Monomial{C} where C
    return Monomial{C}([x, y])
end

function Base.:(*)(x::BasisElement{C}, y::Monomial{C}) :: Monomial{C} where C
    return Monomial{C}([x; unpack(y)])
end

function Base.:(*)(x::Monomial{C}, y::BasisElement{C}) :: Monomial{C} where C
    return Monomial{C}([unpack(x); y])
end

function Base.:(*)(x::Monomial{C}, y::Monomial{C}) :: Monomial{C} where C
    return Monomial{C}([unpack(x); unpack(y)])
end

function Base.:(*)(c::C, b::BasisElement{C}) :: AlgebraElement{C} where C
    return AlgebraElement{C}([(c, Monomial{C}(b))])
end

function Base.:(*)(c::C, m::Monomial{C}) :: AlgebraElement{C} where C
    return AlgebraElement{C}([(c, m)])
    # return AlgebraElement{C}([(c*coeff, mon) for (coeff, mon) in AlgebraElement{C}(m)])
end

function Base.:(*)(c::Union{Int64, C}, a::AlgebraElement{C}) :: AlgebraElement{C} where C
    return AlgebraElement{C}([(c*coeff, mon) for (coeff, mon) in AlgebraElement{C}(a)])
end

function Base.:(*)(c::Int64, b::BasisElement{C}) :: AlgebraElement{C} where C
    return C(c) * b
end

function Base.:(*)(c::Int64, m::Monomial{C}) :: AlgebraElement{C} where C
    return C(c) * m
end

function Base.:(*)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return c * a
end

function Base.:(*)(a1::AlgebraElement{C}, a2::AlgebraElement{C}) :: AlgebraElement{C} where C
    return collectSummands(AlgebraElement{C}([(c1*c2, m1*m2) for (c1, m1) in a1 for (c2, m2) in a2]))
end

function Base.:(*)(a::AlgebraElement{C}, m::Union{BasisElement{C}, Monomial{C}}) :: AlgebraElement{C} where C
    return collectSummands(AlgebraElement{C}([(c, m2*m) for (c, m2) in a]))
end

function Base.:(*)(m::Union{BasisElement{C}, Monomial{C}}, a::AlgebraElement{C}) :: AlgebraElement{C} where C
    return collectSummands(AlgebraElement{C}([(c, m*m2) for (c, m2) in a]))
end


function Base.:(^)(m::Union{BasisElement{C}, Monomial{C}}, n::Int64) :: Monomial{C} where C
    @assert n >= 0
    return prod([m for _ in 1:n]; init=one(Monomial{C}))
end

function Base.:(^)(a::AlgebraElement{C}, n::Int64) :: AlgebraElement{C} where C
    @assert n >= 0
    return prod([a for _ in 1:n])
end


function Base.:(-)(a::StandardOperand{C}) :: AlgebraElement{C} where C
    return (-1)*a
end

function Base.:(-)(x::Operand{C}, y::StandardOperand{C}) :: AlgebraElement{C} where C
    return x + (-y)
end

function Base.:(-)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return a + (-c)
end


comm(x, y) = x*y - y*x


function groupBy(pred, v)
    if isempty(v)
        return typeof(v)()
    end

    tmp = [0; findall([!pred(v[i], v[i+1]) for i in 1:length(v)-1]); length(v)]

    return [v[tmp[i]+1:tmp[i+1]] for i in 1:length(tmp)-1]
end

function prettyPrint(b::BasisElement{C}) :: String where C
    return string(b[1], '(', b[2], ')')
end

function prettyPrint(m::Monomial{C}) :: String where C
    if isone(m)
        return "[]"
    else
        g = groupBy((x,y) -> x[1] === y[1], m)
        innerFormat = map(l -> string(l[1][1], '(', join([string(t[2]) for t in l], ", "), ')'), g)
        return join(innerFormat, "")
    end
end

function prettyPrint(t::Tuple{C, Monomial{C}}) :: String where C
    if isone(t[2])
        return string(isinteger(t[1]) ? Int(t[1]) : t[1])
    elseif isone(t[1])
        return prettyPrint(t[2])
    elseif isone(-t[1])
        return string('-', prettyPrint(t[2]))
    else
        return string(isinteger(t[1]) ? Int(t[1]) : t[1], '⋅', prettyPrint(t[2]))
    end
end

function prettyPrint(a::AlgebraElement{C}) :: String where C
    if iszero(a)
        return string(zero(C))
    else
        return replace(join(map(prettyPrint, a), " + "), "+ -" => "- ")
    end
end

function Base.show(io::IO, ::MIME"text/plain", x::StandardOperand{C}) :: Nothing where C
    print(io, prettyPrint(x))
end
