function symmetrize(f::FreeAssAlgElem)
    R = parent(f)
    g = zero(R)
    for (c, exp) in zip(AbstractAlgebra.coefficients(f), AbstractAlgebra.exponent_words(f))
        g += R(
            fill(c * QQ(1, factorial(length(exp))), factorial(length(exp))),
            [ind for ind in permutations(exp)],
        )

    end
    return g
end
