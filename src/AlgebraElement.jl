include("Structs.jl")

abstract type AbstractAlgebraElement <: Wrapper end

BasisIndex = Int64
BasisElementInternal = Tuple{Symbol, Int64}
MonomialInternal = Vector{BasisElementInternal}
AlgebraElementInternal{C} = Vector{Tuple{C, MonomialInternal}}


struct BasisElement <: AbstractAlgebraElement
    b :: BasisElementInternal

    BasisElement(sym::Symbol, ind::Int64) = new((sym, ind))
    BasisElement(b::BasisElementInternal) = new(b)
    BasisElement(b::BasisElement) = b
end

struct Monomial <: AbstractAlgebraElement
    m :: MonomialInternal

    Monomial() = new([])
    Monomial(b::BasisElementInternal) = new([b])
    Monomial(b::BasisElement) = Monomial(unpack(b))

    Monomial(m::MonomialInternal) = new(m)
    Monomial(m::Vector{BasisElement}) = Monomial(map(unpack, m))
    Monomial(m::Monomial) = m
end

struct AlgebraElement{C} <: AbstractAlgebraElement
    a :: AlgebraElementInternal{C}

    # emtpy list represents empty sum, which is 0
    AlgebraElement{C}() where C = new{C}([])

    AlgebraElement{C}(c::Int64) where C = iszero(c) ? AlgebraElement{C}() : new{C}([(c, Tuple{Symbol, Int64}[])])
    AlgebraElement{C}(c::C) where C = new{C}([(c, Tuple{Symbol, Int64}[])])

    AlgebraElement{C}(b::BasisElementInternal) where C = new{C}([(C(1), [b])])
    AlgebraElement{C}(b::BasisElement) where C = AlgebraElement{C}(unpack(b))

    AlgebraElement{C}(m::MonomialInternal) where C = new{C}([(C(1), m)])
    AlgebraElement{C}(m::Vector{BasisElement}) where C = AlgebraElement{C}(map(unpack, m))
    AlgebraElement{C}(m::Monomial) where C = AlgebraElement{C}(unpack(m))

    AlgebraElement{C}(a::AlgebraElementInternal) where C = new{C}(a)
    AlgebraElement{C}(a::Vector{Tuple{C, Vector{BasisElement}}}) where C = new{C}([(c, map(unpack, m)) for (c, m) in a])
    AlgebraElement{C}(a::Vector{Tuple{C, Monomial}}) where C = new{C}([(c, unpack(m)) for (c, m) in a])
    AlgebraElement{C}(a::AlgebraElement{C}) where C = a
end

StandardOperand{C} = Union{BasisElement, Monomial, AlgebraElement{C},
    BasisElementInternal, MonomialInternal, AlgebraElementInternal{C}}
Operand{C} = Union{Int64, C, BasisElement, Monomial, AlgebraElement{C},
    BasisElementInternal, MonomialInternal, AlgebraElementInternal{C}}


function Base.show(io::IO, b::BasisElement) :: Nothing
    print(io, b[1], "(", b[2], ")")
end

function groupBy(pred, v)
    if isempty(v)
        return typeof(v)()
    end

    tmp = [0; findall([!pred(v[i], v[i+1]) for i in 1:length(v)-1]); length(v)]

    return [v[tmp[i]+1:tmp[i+1]] for i in 1:length(tmp)-1]
end


function Base.show(io::IO, m::Monomial) :: Nothing
    if isempty(m)
        print(io, "[]")
    else
        for g in groupBy((x,y) -> x[1] === y[1], m)
            print(io, g[1][1], "(")
            for i in 1:length(g)-1
                print(io, g[i][2], ", ")
            end
            print(io, g[end][2], ")")
        end
    end
end


function _monomials(a::AlgebraElement) :: Vector{MonomialInternal}
    return unique([mon for (coeff, mon) in a])
end

function monomials(a::AlgebraElement) :: Vector{Monomial}
    return map(Monomial, _monomials(a))
end

function _basisElements(a::AlgebraElement) :: Vector{BasisElementInternal}
    return unique(vcat(_monomials(a)...))
end

function basisElements(a::AlgebraElement) :: Vector{BasisElement}
    return map(BasisElement, _basisElements(a))
end


function _collectSummands(a::Operand{C}) :: AlgebraElementInternal{C} where C
    res = AlgebraElementInternal{C}()

    for (coeff, mon) in AlgebraElement{C}(a)
        i = findfirst(x -> x[2] == mon, res)

        if i === nothing
            push!(res, (coeff, mon))
        else
            res[i] = (res[i][1] + coeff, res[i][2])
        end
    end

    return filter(x -> !iszero(x[1]), res)
end

function collectSummands(a::Operand{C}) :: AlgebraElement{C} where C
    return AlgebraElement{C}(_collectSummands(a))
end

function _collectSummands(as::Vector{AlgebraElement{C}}) :: AlgebraElementInternal{C} where C
    return _collectSummands(vcat(map(unpack, as)...))
end

function collectSummands(as::Vector{AlgebraElement{C}}) :: AlgebraElement{C} where C
    return collectSummands(vcat(map(unpack, as)...))
end


function sameSum(a1::Operand{C}, a2::Operand{C}) :: Bool where C
    return issetequal(_collectSummands(a1), _collectSummands(a2))
end

≐ = sameSum # type symbol via \doteq, autocomplete with tab


function Base.:(+)(x::Operand{C}, as::Vararg{StandardOperand{C}}) :: AlgebraElement{C} where C
    return collectSummands(map(AlgebraElement{C}, [x; as]))
end

function Base.:(+)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return c + a
end


function Base.:(*)(x::Union{BasisElement, Monomial},
                   y::Union{BasisElement, Monomial}) :: Monomial
    return Monomial([unpack(x); unpack(y)])
end

function Base.:(*)(c::Union{Int64, C}, a::StandardOperand{C}) :: AlgebraElement{C} where C
    return AlgebraElement{C}([(c*coeff, mon) for (coeff, mon) in AlgebraElement{C}(a)])
end

function Base.:(*)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return c * a
end

function Base.:(*)(a1::AlgebraElement{C}, a2::AlgebraElement{C}) :: AlgebraElement{C} where C
   return collectSummands([(c1*c2, [m1;m2]) for (c1, m1) in a1 for (c2, m2) in a2])
end

function Base.:(*)(a::AlgebraElement{C}, m::Union{BasisElement, Monomial}) :: AlgebraElement{C} where C
    return a * AlgebraElement{C}(m)
end

function Base.:(*)(m::Union{BasisElement, Monomial}, a::AlgebraElement{C}) :: AlgebraElement{C} where C
    return AlgebraElement{C}(m) * a
end


function Base.:(^)(m::Union{BasisElement, Monomial}, n::Int64) :: Monomial
    @assert n >= 0
    return prod([m for _ in 1:n], init=Monomial())
end

function Base.:(^)(a::AlgebraElement{C}, n::Int64) :: AlgebraElement{C} where C
    @assert n >= 0
    return prod([a for _ in 1:n], init=AlgebraElement{C}(1))
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
