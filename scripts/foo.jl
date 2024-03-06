function foo(L, V, d)
    sp = smash_product(L, V)
    bs = ArcDiagDeformBasis{QQFieldElem}[]
    for i in 0:d
        print(i, "\t")
        push!(bs, ArcDiagDeformBasis{QQFieldElem}(sp, i:i))
        ms = all_pbwdeformations(sp, bs[end])
        print(length(ms), "\t")
        ms = all_pbwdeformations(sp, UnionBasis{QQFieldElem}(bs))
        print(length(ms), "\n")
    end
end
