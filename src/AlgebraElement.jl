BasisIndex = Int64
BasisElement = Tuple{Symbol, Int64}
Monomial{T} = Vector{T}
Scaled{T, C} = Tuple{C, T}
LinearCombination{T, C} = Vector{Scaled{T, C}}
AlgebraElement{C} = LinearCombination{Monomial{BasisElement}, C}
StandardOperand{C} = Union{BasisElement, Monomial{BasisElement}, Scaled{Monomial{BasisElement}, C}, AlgebraElement{C}}
Operand{C} = Union{Int64, C, BasisElement, Monomial{BasisElement}, Scaled{Monomial{BasisElement}, C}, AlgebraElement{C}}


# emtpy list represents empty sum, which is 0
algebraElement(::Type{C}) where C = AlgebraElement{C}()
algebraElement(::Type{C}, c::Union{Int64, C}) where C = iszero(c) ? AlgebraElement() : [(C(c), Monomial{BasisElement}())] :: AlgebraElement{C}
algebraElement(::Type{C}, b::BasisElement) where C = [(C(1), [b])] :: AlgebraElement{C}
algebraElement(::Type{C}, m::Monomial{BasisElement}) where C = [(C(1), m)] :: AlgebraElement{C}
algebraElement(::Type{C}, s::Scaled{Monomial{BasisElement}, C}) where C = [s] :: AlgebraElement{C}
algebraElement(::Type{C}, a::AlgebraElement{C}) where C = a :: AlgebraElement{C}
algebraElement(::Type{C}, x) where C = AlgebraElement{C}(x)
Base.convert(::Type{AlgebraElement{C}}, x::Operand) where C = algebraElement{C}(x)

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
    # castCoeff(c) = isinteger(c) ? Int(c) : c

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

function collectSummands(a::Operand{C}) :: AlgebraElement{C} where C
    res = algebraElement(C)

    for (coeff, mon) in algebraElement(C, a)
        i = findfirst(x -> x[2] == mon, res)

        if i === nothing
            push!(res, (coeff, mon))
        else
            res[i] = (res[i][1] + coeff, res[i][2])
        end
    end

    return filter(x -> !iszero(x[1]), res)
end

function collectSummands(as::Vector{AlgebraElement{C}}) :: AlgebraElement{C} where C
    return collectSummands(algebraElement(C, vcat(as...)))
end

function sameSum(a1::Operand{C}, a2::Operand{C}) :: Bool where C
    return issetequal(collectSummands(a1), collectSummands(a2))
end

≐ = sameSum # type symbol via \doteq, autocomplete with tab


function Base.:(+)(x::Operand{C}, as::Vararg{StandardOperand{C}}) :: AlgebraElement{C} where C
    return collectSummands(map(y -> algebraElement(C, y), [x, as...]))
end

function Base.:(+)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return c + a
end


function Base.:(*)(x::Union{BasisElement, Monomial{BasisElement}},
                   y::Union{BasisElement, Monomial{BasisElement}}) :: Monomial{BasisElement}
    return [x; y]
end

function Base.:(*)(c::Union{Int64, C}, a::StandardOperand{C}) :: AlgebraElement{C} where C
    return [(c*coeff, mon) for (coeff, mon) in algebraElement(C, a)]
end

function Base.:(*)(a::StandardOperand{C}, c::Union{Int64, C}) :: AlgebraElement{C} where C
    return c * a
end

function Base.:(*)(a1::Union{Scaled{Monomial{BasisElement}, C}, AlgebraElement{C}},
                   a2::Union{Scaled{Monomial{BasisElement}, C}, AlgebraElement{C}}) :: AlgebraElement{C} where C
   return collectSummands([(c1*c2, [m1;m2]) for (c1, m1) in algebraElement(C, a1) for (c2, m2) in algebraElement(C, a2)])
end

function Base.:(*)(a::Union{Scaled{Monomial{BasisElement}, C}, AlgebraElement{C}},
                   m::Union{BasisElement, Monomial{BasisElement}}) :: AlgebraElement{C} where C
    return a * algebraElement(C, m)
end

function Base.:(*)(m::Union{BasisElement, Monomial{BasisElement}}, a::AlgebraElement{C}) :: AlgebraElement{C} where C
    return algebraElement(C, m) * a
end


function Base.:(^)(m::Union{BasisElement, Monomial{BasisElement}}, n::Int64) :: Monomial{BasisElement} where C
    @assert n >= 0
    return prod([m for _ in 1:n], init=Monomial{BasisElement}())
end

function Base.:(^)(a::Union{Scaled{Monomial{BasisElement}}, AlgebraElement{C}}, n::Int64) :: AlgebraElement{C} where C
    @assert n >= 0
    return prod([a for _ in 1:n], init=algebraElement(C, 1))
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
