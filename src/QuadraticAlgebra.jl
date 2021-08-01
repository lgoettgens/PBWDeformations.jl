struct QuadraticAlgebra{C, T}
    basis :: Vector{BasisElement{C}}

    """
    Stores relations of the form ab = c for basis elements a,b.
    An empty list represents the empty sum, in which case ab = 0.
    An absent entry means that there is no relation, so we cannot simplify ab.
    """
    relTable :: Dict{Tuple{BasisElement{C}, BasisElement{C}}, AlgebraElement{C}}
    extraData :: T

    QuadraticAlgebra{C, T}(basis, relTable, extraData = nothing) where {C, T} =
        new{C, T}(basis, relTable, extraData)
end

function Base.:(==)(alg1::QuadraticAlgebra{C}, alg2::QuadraticAlgebra{C}) :: Bool where C
    (alg1.basis, alg1.relTable, alg1.extraData) ==
    (alg2.basis, alg2.relTable, alg2.extraData)
end

function Base.show(io::IO, alg::QuadraticAlgebra{C}) :: Nothing where C
    println(io, "Algebra with quadratic relations of dimension ", length(alg.basis))
    println(io, "Relation table has ", length(alg.relTable), " entries")
    println(io, "Coefficients of type ", C)
    println(io, "Extra data:")
    print(io, alg.extraData)
end

function Base.in(b::BasisElement{C}, alg::QuadraticAlgebra{C}) :: Bool where C
    return b in alg.basis
end

function Base.in(m::Monomial{C}, alg::QuadraticAlgebra{C}) :: Bool where C
    return all(b -> b in alg, m)
end

function Base.in(a::AlgebraElement{C}, alg::QuadraticAlgebra{C}) :: Bool where C
    return all(m in alg for (c, m) in a)
end


function normalForm(alg::QuadraticAlgebra{C}, a::AlgebraElement{C}) :: AlgebraElement{C} where C
    todo = copy(unpack(a))
    result = AlgebraElement{C}(0)

    while !isempty(todo)
        coeff, mon = pop!(todo)

        changed = false
        for i in 1:length(mon)-1
            if haskey(alg.relTable, (mon[i], mon[i+1]))
                changed = true

                # TODO: something like this: todo += coeff * (mon[1:i-1] * alg.relTable[(mon[i], mon[i+1])] * mon[i+2:end])
                todo = unpack(AlgebraElement{C}(todo) + coeff * (Monomial{C}(mon[1:i-1]) * alg.relTable[(mon[i], mon[i+1])] * Monomial{C}(mon[i+2:end])))

                break
            end
        end

        if !changed
            result += coeff * mon
        end
    end

    return result
end

function normalForm(alg::QuadraticAlgebra{C}, m::Union{BasisElement{C}, Monomial{C}}) :: AlgebraElement{C} where C
    return normalForm(alg, AlgebraElement{C}(m))
end
