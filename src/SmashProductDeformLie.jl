mutable struct SmashProductDeformLie{C <: RingElement}
    dimL :: Int64
    dimV :: Int64
    baseL :: Vector{QuadraticQuoAlgebraElem{C}}
    baseV :: Vector{QuadraticQuoAlgebraElem{C}}
    coeff_ring :: Ring
    alg :: QuadraticQuoAlgebra{C}
    # dynkin :: Char
    # n :: Int64
    # lambda :: Vector{Int64}
    # matrixRepL :: Vector{Matrix{Int64}}
    symmetric :: Bool
    kappa :: Matrix{QuadraticQuoAlgebraElem{C}}
end


function smash_product_deform_lie(sp::SmashProductLie{C}, kappa::Matrix{QuadraticQuoAlgebraElem{C}}) where C <: RingElement
    size(kappa) == (sp.dimV, sp.dimV) || throw(ArgumentError("kappa has wrong dimensions."))
    
    dimL = sp.dimL
    dimV = sp.dimV
    coeff_ring = sp.coeff_ring
    baseL = sp.baseL
    baseV = sp.baseV

    for i in 1:dimV, j in 1:i
        kappa[i,j] == -kappa[j,i] || throw(ArgumentError("kappa is not skew-symmetric."))
        all(x -> x <= dimL, var_ids(kappa[i,j])) || throw(ArgumentError("kappa does not only take values in the hopf algebra"))
        all(x -> x <= dimL, var_ids(kappa[j,i])) || throw(ArgumentError("kappa does not only take values in the hopf algebra"))
    end

    symmetric = true
    rels = Dict{Tuple{Int,Int}, QuadraticQuoAlgebraElem{C}}()
    for i in 1:dimV, j in 1:dimV
        # We have the commutator relation [v_i, v_j] = kappa[i,j]
        # which is equivalent to v_i*v_j = v_j*v_i + kappa[i,j]
        rels[(dimL+i, dimL+j)] = baseV[j] * baseV[i] + kappa[i,j]
        symmetric &= iszero(kappa[i,j])
    end

    alg = quadratic_quo_algebra(sp.alg, rels)
    baseL = map(alg, baseL)
    baseV = map(alg, baseV)

    return SmashProductDeformLie{C}(dimL, dimV, baseL, baseV, coeff_ring, alg, symmetric, kappa), (baseL, baseV)
end

function smash_product_symmdeform_lie(sp::SmashProductLie{C}) where C <: RingElement
    kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
    return smash_product_deform_lie(sp, kappa)
end

ngens(d::SmashProductDeformLie) = d.dimL, d.dimV

function gens(d::SmashProductDeformLie{C}) where C <: RingElement
    return [gen(d.alg, i) for i in 1:d.dimL], [gen(d.alg, i+d.dimL) for i in 1:d.dimV]
end


function show(io::IO, deform::SmashProductDeformLie)
    local max_gens = 4 # largest number of generators to print
    if deform.symmetric
        print(io, "Symmetric ")
    end
    print(io, "Deformation of ")
    print(io, "Lie Algebra Smash Product with basis ")
    for i = 1:min(deform.dimL - 1, max_gens - 1)
       print(io, string(deform.alg.S[i]), ", ")
    end
    if deform.dimL > max_gens
       print(io, "..., ")
    end
    print(io, string(deform.alg.S[deform.dimL]) * ", ")
    for i = 1:min(deform.dimV - 1, max_gens - 1)
        print(io, string(deform.alg.S[deform.dimL+i]), ", ")
     end
     if deform.dimV > max_gens
        print(io, "..., ")
     end
     print(io, string(deform.alg.S[deform.dimL+deform.dimV]))
    print(io, " over ")
    print(IOContext(io, :compact => true), deform.coeff_ring)
end

function change_base_ring(R::Ring, d::SmashProductDeformLie{C}) where C <: RingElement
    alg = change_base_ring(R, d.alg)
    baseL = [gen(alg, i) for i in 1:d.dimL]
    baseV = [gen(alg, d.dimL+i) for i in 1:d.dimV]
    kappa = map(alg, d.kappa)

    return SmashProductDeformLie{elem_type(R)}(d.dimL, d.dimV, baseL, baseV, R, alg, d.symmetric, kappa)
end


function pbwdeform_eqs(deform::SmashProductDeformLie{C}) where C <: RingElement
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras.
    # DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011
    dimL = deform.dimL
    dimV = deform.dimV
    kappa = deform.kappa
    x(i) = gen(deform.alg, i)
    v(i) = gen(deform.alg, dimL + i)

    @debug "Equation generation, first phase"

    ## (a) κ is H-invariant
    iter_a = (comm(comm(h, v(i), true), v(j)) # κ([h⋅v_i,v_j])
            + comm(v(i), comm(h, v(j), true)) # κ([v_i,h⋅v_j])
            - comm(h, comm(v(i), v(j), true))  # h⋅κ([v_i,v_j])
        for (h, (i,j)) in Iterators.product([x(i) for i in 1:dimL], Combinatorics.Combinations(dimV, 2)))
    # m[1][2] denotes the index of the only basis element in the monomial m

    ## (b) trivial
    iter_b = []

    @debug "Equation generation, second phase"

    ## (c) 0 = κ ⊗ id - id ⊗ κ on (I ⊗ V) ∩ (V ⊗ I)
    # (I ⊗ V) ∩ (V ⊗ I) has basis v_iv_jv_k + v_jv_kv_i + v_kv_iv_j - v_kv_jv_i - v_jv_iv_k - v_iv_kv_j for i<j<k
    iter_c = (comm(v(i), v(j), true)*v(k) - v(i)*comm(v(j), v(k), true)
            + comm(v(j), v(k), true)*v(i) - v(j)*comm(v(k), v(i), true)
            + comm(v(k), v(i), true)*v(j) - v(k)*comm(v(i), v(j), true)
            - comm(v(k), v(j), true)*v(i) + v(k)*comm(v(j), v(i), true)
            - comm(v(j), v(i), true)*v(k) + v(j)*comm(v(i), v(k), true)
            - comm(v(i), v(k), true)*v(j) + v(i)*comm(v(k), v(j), true)
        for (i,j,k) in Combinatorics.Combinations(dimV, 3))

    ## (d) trivial
    iter_d = []

    iter = Iterators.flatten([iter_a, iter_b, iter_c, iter_d])
    return Iterators.map(normal_form, iter)
end

function pbwdeform_neqs(deform::SmashProductDeformLie{C}) where C <: RingElement
    num_a = deform.dimL * binomial(deform.dimV, 2)
    num_b = 0
    num_c = binomial(deform.dimV, 3)
    num_d = 0

    return num_a + num_b + num_c + num_d
end

function ispbwdeform(d::SmashProductDeformLie{C}) :: Bool where C <: RingElement
    return all(iszero, pbwdeform_eqs(d))
end
