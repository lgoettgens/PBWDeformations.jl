BasisIndex = Int64
Coefficient = Rational{Int64}
BasisElement = Tuple{Symbol, Int64}
Monomial{T} = Vector{T}
Scaled{T} = Tuple{Coefficient, T}
LinearCombination{T} = Vector{Scaled{T}}
AlgebraElement = LinearCombination{Monomial{BasisElement}}
StandardOperand = Union{BasisElement, Monomial{BasisElement}, Scaled{Monomial{BasisElement}}, AlgebraElement}
Operand = Union{Int64, Coefficient, BasisElement, Monomial{BasisElement}, Scaled{Monomial{BasisElement}}, AlgebraElement}


# emtpy list represents empty sum, which is 0
algebraElement() = AlgebraElement()
algebraElement(c::Union{Int64, Coefficient}) = iszero(c) ? AlgebraElement() : [(Coefficient(c), Monomial{BasisElement}())] :: AlgebraElement
algebraElement(b::BasisElement) = [(Coefficient(1), [b])] :: AlgebraElement
algebraElement(m::Monomial{BasisElement}) = [(Coefficient(1), m)] :: AlgebraElement
algebraElement(s::Scaled{Monomial{BasisElement}}) = [s] :: AlgebraElement
algebraElement(a::AlgebraElement) = a :: AlgebraElement
algebraElement(x) = AlgebraElement(x)
Base.convert(::Type{AlgebraElement}, x::Operand) = algebraElement(x)

function groupBy(pred, v::Vector{T}) where T
    if isempty(v)
        return Vector{T}[]
    end

    tmp = [0; findall([!pred(v[i], v[i+1]) for i in 1:length(v)-1]); length(v)]

    return [v[tmp[i]+1:tmp[i+1]] for i in 1:length(tmp)-1]
end

function Base.show(io::IO, b::BasisElement) :: Nothing
    print(io, b[1], "(", b[2], ")")
end

function Base.show(io::IO, m::Monomial{BasisElement}) :: Nothing
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

function Base.show(io::IO, s::Scaled{Monomial{BasisElement}}) :: Nothing
    castCoeff(c) = isinteger(c) ? Int(c) : c

    if iszero(s[1])
        print(io, 0)
    elseif isempty(s[2])
        print(io, castCoeff(s[1]))
    elseif isone(s[1])
        print(io, s[2])
    else
        print(io, castCoeff(s[1]), "⋅", s[2])
    end
end


function monomials(a::AlgebraElement) :: Vector{Monomial{BasisElement}}
    return unique([mon for (coeff, mon) in a])
end

function basisElements(a::AlgebraElement) :: Vector{BasisElement}
    return unique(vcat(monomials(a)...))
end

function collectSummands(a::Operand) :: AlgebraElement
    res = algebraElement()

    for (coeff, mon) in algebraElement(a)
        i = findfirst(x -> x[2] == mon, res)

        if i === nothing
            push!(res, (coeff, mon))
        else
            res[i] = (res[i][1] + coeff, res[i][2])
        end
    end

    return filter(x -> !iszero(x[1]), res)
end

function collectSummands(as::Vector{AlgebraElement}) :: AlgebraElement
    return collectSummands(algebraElement(vcat(as...)))
end

function sameSum(a1::Operand, a2::Operand) :: Bool
    return issetequal(collectSummands(a1), collectSummands(a2))
end

≐ = sameSum # type symbol via \doteq, autocomplete with tab


function Base.:(+)(x::Operand, as::Vararg{StandardOperand}) :: AlgebraElement
    return collectSummands(map(algebraElement, [x, as...]))
end

function Base.:(+)(a::StandardOperand, c::Union{Int64, Coefficient}) :: AlgebraElement
    return c + a
end


function Base.:(*)(x::Union{BasisElement, Monomial{BasisElement}},
                   y::Union{BasisElement, Monomial{BasisElement}}) :: Monomial{BasisElement}
    return [x; y]
end

function Base.:(*)(c::Union{Int64, Coefficient}, a::StandardOperand) :: AlgebraElement
    return [(c*coeff, mon) for (coeff, mon) in algebraElement(a)]
end

function Base.:(*)(a::StandardOperand, c::Union{Int64, Coefficient}) :: AlgebraElement
    return c * a
end

function Base.:(*)(a1::Union{Scaled{Monomial{BasisElement}}, AlgebraElement},
                   a2::Union{Scaled{Monomial{BasisElement}}, AlgebraElement}) :: AlgebraElement
   return collectSummands([(c1*c2, [m1;m2]) for (c1, m1) in algebraElement(a1) for (c2, m2) in algebraElement(a2)])
end

function Base.:(*)(a::Union{Scaled{Monomial{BasisElement}}, AlgebraElement},
                   m::Union{BasisElement, Monomial{BasisElement}}) :: AlgebraElement
    return a * algebraElement(m)
end

function Base.:(*)(m::Union{BasisElement, Monomial{BasisElement}}, a::AlgebraElement) :: AlgebraElement
    return algebraElement(m) * a
end


function Base.:(^)(m::Union{BasisElement, Monomial{BasisElement}}, n::Int64) :: Monomial{BasisElement}
    @assert n >= 0
    return prod([m for _ in 1:n], init=Monomial{BasisElement}())
end

function Base.:(^)(a::Union{Scaled{Monomial{BasisElement}}, AlgebraElement}, n::Int64) :: AlgebraElement
    @assert n >= 0
    return prod([a for _ in 1:n], init=algebraElement(1))
end


function Base.:(-)(a::StandardOperand) :: AlgebraElement
    return (-1)*a
end

function Base.:(-)(x::Operand, y::StandardOperand) :: AlgebraElement
    return x + (-y)
end

function Base.:(-)(a::StandardOperand, c::Union{Int64, Coefficient}) :: AlgebraElement
    return a + (-c)s
end


comm(x, y) = x*y - y*x
