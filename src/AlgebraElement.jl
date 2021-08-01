abstract type AbstractAlgebraElement <: Wrapper end

BasisIndex = Int64

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
end

struct AlgebraElement{C} <: AbstractAlgebraElement
    a :: Vector{Tuple{C, Monomial{C}}}

    # emtpy list represents empty sum, which is 0
    AlgebraElement{C}() where C = new{C}([])

    AlgebraElement{C}(c::Int64) where C = iszero(c) ? AlgebraElement{C}() : new{C}([(C(c), Monomial{C}())])
    AlgebraElement{C}(c::C) where C = iszero(c) ? AlgebraElement{C}() : new{C}([(c, Monomial{C}())])

    AlgebraElement{C}(b::BasisElement{C}) where C = new{C}([(C(1), Monomial{C}(b))])

    AlgebraElement{C}(m::Monomial{C}) where C = new{C}([(C(1), m)])

    AlgebraElement{C}(a::Vector{Tuple{C, Monomial{C}}}) where C = new{C}(a)
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
Operand{C} = Union{Int64, C, BasisElement{C}, Monomial{C}, AlgebraElement{C}}

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


function monomials(a::AlgebraElement{C}) :: Vector{Monomial{C}} where C
    return unique([mon for (coeff, mon) in a])
end

function basisElements(a::AlgebraElement{C}) :: Vector{BasisElement{C}} where C
    # TODO: implement vcat functionality for wrapper structs
    return unique(vcat(map(unpack, monomials(a))...))
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
    # TODO: implement vcat functionality for wrapper structs
    return collectSummands(AlgebraElement{C}(vcat(map(unpack, as)...)))
end

function sameSum(a1::Operand{C}, a2::Operand{C}) :: Bool where C
    # https://github.com/JuliaLang/julia/issues/41748 
    _issubset(a, b) = all(x -> x in b, a) # for vectors of length >70 issubset is weird
    _issetequal(a, b) = _issubset(a, b) && _issubset(b, a)

    return _issetequal(collectSummands(AlgebraElement{C}(a1)),
                      collectSummands(AlgebraElement{C}(a2)))
end

≐ = sameSum # type symbol via \doteq, autocomplete with tab

function Base.iszero(a::AlgebraElement{C}) :: Bool where C
    return a ≐ 0
end


function Base.:(+)(x::Operand{C}, y::StandardOperand{C}) :: AlgebraElement{C} where C
    return collectSummands([AlgebraElement{C}(x), AlgebraElement{C}(y)])
end

function Base.:(+)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return c + a
end


function Base.:(*)(x::BasisElement{C}, y::BasisElement{C}) :: Monomial{C} where C
    return Monomial{C}([x; y])
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

function Base.:(*)(c::Union{Int64, C}, a::StandardOperand{C}) :: AlgebraElement{C} where C
    return AlgebraElement{C}([(c*coeff, mon) for (coeff, mon) in AlgebraElement{C}(a)])
end

function Base.:(*)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return c * a
end

function Base.:(*)(a1::AlgebraElement{C}, a2::AlgebraElement{C}) :: AlgebraElement{C} where C
    # TODO ; for wrapper 
   return collectSummands(AlgebraElement{C}([(c1*c2, m1*m2) for (c1, m1) in a1 for (c2, m2) in a2]))
end

function Base.:(*)(a::AlgebraElement{C}, m::Union{BasisElement{C}, Monomial{C}}) :: AlgebraElement{C} where C
    return a * AlgebraElement{C}(m)
end

function Base.:(*)(m::Union{BasisElement{C}, Monomial{C}}, a::AlgebraElement{C}) :: AlgebraElement{C} where C
    return AlgebraElement{C}(m) * a
end


function Base.:(^)(m::Union{BasisElement{C}, Monomial{C}}, n::Int64) :: Monomial{C} where C
    @assert n >= 0
    return prod([m for _ in 1:n]; init=Monomial{C}())
end

function Base.:(^)(a::AlgebraElement{C}, n::Int64) :: AlgebraElement{C} where C
    @assert n >= 0
    return prod([a for _ in 1:n]; init=AlgebraElement{C}(1))
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
