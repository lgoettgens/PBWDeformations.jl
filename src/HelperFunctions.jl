BasisIndex = Int64
Coefficient = Rational
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

Base.convert(::Type{Product{BasisElement}}, b::BasisElement) = [b]
Base.convert(::Type{AlgebraElement}, b::BasisElement) = [(1, [b])]
Base.convert(::Type{AlgebraElement}, p::Product{BasisElement}) = [(1, p)]
Base.convert(::Type{AlgebraElement}, c::Coefficient) = [(c, [])]

function Base.:(+)(a1::AlgebraElement, a2::AlgebraElement) :: AlgebraElement
    res = copy(a1)

    for (coeff, prod) in a2
        i = findfirst(x -> x[2] == prod, res)

        if i === nothing
            append!(res, (coeff, prod))
        else
            res[i][1] += coeff
        end
    end

    return res
end

function Base.:(+)(a1::Union{Coefficient, BasisElement, Product{BasisElement}, AlgebraElement},
                   a2::Union{Coefficient, BasisElement, Product{BasisElement}, AlgebraElement}) :: AlgebraElement
    return convert(AlgebraElement, a1) + convert(AlgebraElement, a2)
end

function Base.:(+)(l1::LinearCombination{Vector{BasisIndex}}, l2::LinearCombination{Vector{BasisIndex}}) :: LinearCombination{Vector{BasisIndex}}
    res = copy(l1)

    for (coeff, ind) in l2
        i = findfirst(x -> x[2] == ind, res)

        if i === nothing
            append!(res, (coeff, ind))
        else
            res[i][1] += coeff
        end
    end

    return res
end


function Base.:(*)(p::Product{BasisElement}, x::Union{BasisElement, Product{BasisElement}}) :: Product{BasisElement}
    return append!(p, x)
end

function Base.:(*)(x::Union{BasisElement, Product{BasisElement}}, p::Product{BasisElement}) :: Product{BasisElement}
    return prepend!(p, x)
end

function Base.:(*)(a1::AlgebraElement, a2::AlgebraElement) :: AlgebraElement
    return [(c1*c2, p1*p2) for (c1, p1) in a1 for (c2, p2) in a2]
end

function Base.:(*)(x::Union{Coefficient, BasisElement, Product{BasisElement}}, a::AlgebraElement) :: AlgebraElement
    return convert(AlgebraElement, x) * a
end

function Base.:(*)(a::AlgebraElement, x::Union{Coefficient, BasisElement, Product{BasisElement}}) :: AlgebraElement
    return a * convert(AlgebraElement, x)
end

function Base.:(*)(coeff::Coefficient, x::Union{BasisElement, Product{BasisElement}}) :: AlgebraElement
    return coeff * convert(AlgebraElement, x)
end

function Base.:(*)(x::Union{BasisElement, Product{BasisElement}}, coeff::Coefficient) :: AlgebraElement
    return coeff * x
end
