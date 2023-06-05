

function change_base_ring(R::Ring, sp::SmashProductLie{C, CL}) where {C <: RingElem, CL <: RingElem}
    alg, _ = free_associative_algebra(R, symbols(sp.alg))
    rels = QuadraticRelations{elem_type(R)}(k => change_base_ring(R, a, parent=alg) for (k, a) in sp.rels)

    return SmashProductLie{elem_type(R), CL}(R, sp.L, sp.V, alg, rels)
end

function change_base_ring(R::Ring, d::SmashProductLieDeform{C, CL}) where {C <: RingElem, CL <: RingElem}
    sp = change_base_ring(R, d.sp)
    rels = QuadraticRelations{elem_type(R)}(k => change_base_ring(R, a, parent=sp.alg) for (k, a) in d.rels)
    kappa = map(e -> change_base_ring(R, e, parent=sp.alg), d.kappa)

    return SmashProductLieDeform{elem_type(R), CL}(sp, rels, kappa)
end
