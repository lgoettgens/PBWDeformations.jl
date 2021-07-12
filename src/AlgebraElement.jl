BasisIndex = Int64
Coefficient = Rational{Int64}
BasisElement = Tuple{Symbol, Int64}
Product{T} = Vector{T}
LinearCombination{T} = Vector{Tuple{Coefficient, T}}
AlgebraElement = LinearCombination{Product{BasisElement}}

function monomials(a::AlgebraElement) :: Vector{Product{BasisElement}}
    return unique([prod for (coeff, prod) in a])
end

function basisElements(a::AlgebraElement) :: Vector{BasisElement}
    return unique(vcat(monomials(a)...))
end

function collectSummands(a::Union{AlgebraElement, Vector{AlgebraElement}}) :: AlgebraElement
    a = isa(a, AlgebraElement) ? a : AlgebraElement(vcat(a...))
    res = AlgebraElement()

    for (coeff, prod) in a
        i = findfirst(x -> x[2] == prod, res)

        if i === nothing
            push!(res, (coeff, prod))
        else
            # res[i][1] += coeff does not work since tuples are immutable
            res[i] = (res[i][1] + coeff, res[i][2])
        end
    end

    return filter(x -> !iszero(x[1]), res)
end

function sameSum(a1::AlgebraElement, a2::AlgebraElement) :: Bool
    return issetequal(collectSummands(a1), collectSummands(a2))
end


# emtpy list represents empty sum, which is 0
algebraElement() = AlgebraElement()
algebraElement(c::Union{Int64, Coefficient}) = iszero(c) ? AlgebraElement() : [(Coefficient(c), Product{BasisElement}())] :: AlgebraElement
algebraElement(b::BasisElement) = [(Coefficient(1), [b])] :: AlgebraElement
algebraElement(p::Product{BasisElement}) = [(Coefficient(1), p)] :: AlgebraElement
algebraElement(a::AlgebraElement) = a :: AlgebraElement


function Base.:(+)(x::Union{Int64, Coefficient, BasisElement, Product{BasisElement}, AlgebraElement},
                   as::Vararg{Union{BasisElement, Product{BasisElement}, AlgebraElement}}) :: AlgebraElement
    return collectSummands(map(algebraElement, [x, as...]))
end

function Base.:(+)(a::Union{BasisElement, Product{BasisElement}, AlgebraElement}, c::Union{Int64, Coefficient}) :: AlgebraElement
    return c + a
end


function Base.:(*)(x::Union{BasisElement, Product{BasisElement}},
                   y::Union{BasisElement, Product{BasisElement}}) :: Product{BasisElement}
    return [x; y]
end

function Base.:(*)(c::Union{Int64, Coefficient}, a::Union{BasisElement, Product{BasisElement}, AlgebraElement}) :: AlgebraElement
    return [(c*coeff, mon) for (coeff, mon) in algebraElement(a)]
end

function Base.:(*)(a::Union{BasisElement, Product{BasisElement}, AlgebraElement}, c::Union{Int64, Coefficient}) :: AlgebraElement
    return c * a
end

function Base.:(*)(a1::AlgebraElement, a2::AlgebraElement) :: AlgebraElement
   return collectSummands([(c1*c2, [m1;m2]) for (c1, m1) in a1 for (c2, m2) in a2])
end

function Base.:(*)(a::AlgebraElement, m::Union{BasisElement, Product{BasisElement}}) ::AlgebraElement
    return a * algebraElement(m)
end

function Base.:(*)(m::Union{BasisElement, Product{BasisElement}}, a::AlgebraElement) ::AlgebraElement
    return algebraElement(m) * a
end


function Base.:(^)(m::Union{BasisElement, Product{BasisElement}}, n::Int64) :: Product{BasisElement}
    @assert n >= 0
    return prod([m for _ in 1:n])
end

function Base.:(^)(a::AlgebraElement, n::Int64) :: AlgebraElement
    @assert n >= 0
    return prod([a for _ in 1:n], init=algebraElement(1))
end


function Base.:(-)(a::Union{BasisElement, Product{BasisElement}, AlgebraElement}) :: AlgebraElement
    return (-1)*a
end

function Base.:(-)(x::Union{Int64, Coefficient, BasisElement, Product{BasisElement}, AlgebraElement},
                   y::Union{BasisElement, Product{BasisElement}, AlgebraElement}) :: AlgebraElement
    return x + (-y)
end

function Base.:(-)(a::Union{BasisElement, Product{BasisElement}, AlgebraElement}, c::Union{Int64, Coefficient}) :: AlgebraElement
    return a + (-c)
end
