"""
    pbwdeform_eqs(d::SmashProductLieDeform; disabled::Vector{Symbol}=Symbol[])

Returns the equations for `d` being a PBW deformation of a smash product as
in Theorem 3.1 of [WW14](@cite).
Subsets of the equations can be disabled by passing the corresponding symbols as
keyword arguments, e.g. `disabled = [:c, :d]`.
"""
function pbwdeform_eqs(d::SmashProductLieDeform; disabled::Vector{Symbol}=Symbol[])
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras.
    # DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011
    dimL = dim(lie_algebra(d))
    dimV = dim(lie_module(d))

    x(i) = gen(d, i, :L)
    v(i) = gen(d, i, :V)

    scomm(a, b) = simplify(comm(a, b))

    ## (a) κ is H-invariant
    iter_a =
        :a in disabled ? [] :
        (
            comm(scomm(h, v(i)), v(j))   # κ([h⋅v_i,v_j])
            + comm(v(i), scomm(h, v(j))) # κ([v_i,h⋅v_j])
            - comm(h, scomm(v(i), v(j))) # h⋅κ([v_i,v_j])
            for (h, (i, j)) in Iterators.product([x(i) for i in 1:dimL], combinations(dimV, 2))
        )

    ## (b) trivial
    iter_b = :b in disabled ? [] : []

    ## (c) 0 = κ ⊗ id - id ⊗ κ on (I ⊗ V) ∩ (V ⊗ I)
    # (I ⊗ V) ∩ (V ⊗ I) has basis v_iv_jv_k + v_jv_kv_i + v_kv_iv_j - v_kv_jv_i - v_jv_iv_k - v_iv_kv_j for i<j<k
    iter_c =
        :c in disabled ? [] :
        (
            (scomm(v(i), v(j)) * v(k) - v(i) * scomm(v(j), v(k))) +
            (scomm(v(j), v(k)) * v(i) - v(j) * scomm(v(k), v(i))) +
            (scomm(v(k), v(i)) * v(j) - v(k) * scomm(v(i), v(j))) -
            (scomm(v(k), v(j)) * v(i) - v(k) * scomm(v(j), v(i))) -
            (scomm(v(j), v(i)) * v(k) - v(j) * scomm(v(i), v(k))) -
            (scomm(v(i), v(k)) * v(j) - v(i) * scomm(v(k), v(j))) for (i, j, k) in combinations(dimV, 3)
        )

    ## (d) trivial
    iter_d = :d in disabled ? [] : []

    return Iterators.flatten([iter_a, iter_b, iter_c, iter_d])
end

function pbwdeform_neqs(d::SmashProductLieDeform)
    dimL = dim(lie_algebra(d))
    dimV = dim(lie_module(d))

    num_a = dimL * binomial(dimV, 2)
    num_b = 0
    num_c = binomial(dimV, 3)
    num_d = 0

    return num_a + num_b + num_c + num_d
end

"""
    is_pbwdeformation(d::SmashProductLieDeform)

Check if `d` is a Poincare-Birkhoff-Witt deformation of a smash product.
Uses Theorem 3.1 of [WW14](@cite).
"""
function is_pbwdeformation(d::SmashProductLieDeform)
    return all(iszero, pbwdeform_eqs(d))
end


###############################################################################
#
#   All PBW deformations
#
###############################################################################

function coefficient_comparison(eq::FreeAssAlgElem)
    return collect(coefficients(eq))
end

function linpoly_to_spvector(a::QQMPolyRingElem)
    @assert total_degree(a) <= 1

    return sparse_row(QQ, [(findfirst(==(1), e), c) for (e, c) in zip(exponents(a), coefficients(a))])
end


"""
    all_pbwdeformations(sp::SmashProductLie{C}, deform_basis::DeformBasis{C}; special_return=Nothing) where {C <: RingElem}

Computes a basis of all Poincare-Birkhoff-Witt deformations of `sp`.
`deform_basis` specifies the basis to use for the space of deformation maps.
If `special_return` is `SMat`, the function returns intermediate results.

Uses [`pbwdeform_eqs`](@ref) and thus Theorem 3.1 of [WW14](@cite).
"""
function all_pbwdeformations(
    sp::SmashProductLie{C},
    deform_basis::DeformBasis{C};
    special_return::Type{T}=Nothing,
) where {C <: RingElem, T <: Union{Nothing, SMat}}
    @req coefficient_ring(sp) == QQ "Only implemented for QQ coefficients."

    dimL = dim(lie_algebra(sp))
    dimV = dim(lie_module(sp))

    nvars = length(deform_basis)

    @vprintln :PBWDeformations 1 "Constructing MPolyRing..."
    R, vars = polynomial_ring(coefficient_ring(sp), max(nvars, 1))

    @vprintln :PBWDeformations 1 "Changing SmashProductLie coeffcient type..."
    new_sp = smash_product(R, lie_algebra(sp), lie_module(sp))

    @vprintln :PBWDeformations 1 "Constructing kappa..."
    kappa = fill(zero(underlying_algebra(new_sp)), dimV, dimV)
    for (i, b) in enumerate(deform_basis)
        kappa += vars[i] .* map(e -> change_base_ring(R, e, parent=underlying_algebra(new_sp)), b)
    end

    @vprintln :PBWDeformations 1 "Constructing deformation..."
    d = deform(new_sp, kappa)

    @vprintln :PBWDeformations 1 "Generating equation iterator..."
    neqs = pbwdeform_neqs(d)
    iter = Iterators.map(
        linpoly_to_spvector,
        Iterators.flatten(
            Iterators.map(
                function (x)
                    i = x[1]
                    a = x[2]
                    @vprintln :PBWDeformations 2 "Equations $(lpad(floor(Int, 100*i / neqs), 3))%, $(lpad(i, ndigits(neqs)))/$(neqs)"
                    coefficient_comparison(simplify(a).alg_elem)
                end,
                enumerate(pbwdeform_eqs(d)),
            ),
        ),
    )

    @vprintln :PBWDeformations 1 "Computing reduced row-echelon form..."
    lgs = sparse_matrix(coefficient_ring(sp), 0, nvars)
    for v in iter
        Hecke._add_row_to_rref!(lgs, v)
    end

    @vprintln :PBWDeformations 1 "Computing the kernel..."
    kernel_dim, kernel = right_kernel(lgs)

    if special_return <: SMat
        return kernel, vars
    end

    @vprintln :PBWDeformations 1 "Computing a basis..."
    kappas = Vector{DeformationMap{C}}(undef, kernel_dim)
    for l in 1:kernel_dim
        kappa = fill(zero(underlying_algebra(sp)), dimV, dimV)
        for (i, b) in enumerate(deform_basis)
            kappa += kernel[i, l] .* b
        end
        kappas[l] = kappa
    end
    return kappas
end

"""
    all_pbwdeformations(sp::SmashProductLie{C}, degs::AbstractVector{Int}, DeformBasisType::Type{<:DeformBasis{C}}=StdDeformBasis{C}; special_return=Nothing) where {C <: RingElem}

Computes a basis of all Poincare-Birkhoff-Witt deformations of `sp` of degrees `degs`.
`DeformBasisType` specifies the type of basis to use for the space of deformation maps.
If `special_return` is `SMat`, the function returns intermediate results.

Uses [`pbwdeform_eqs`](@ref) and thus Theorem 3.1 of [WW14](@cite).
"""
function all_pbwdeformations(
    sp::SmashProductLie{C},
    degs::AbstractVector{Int},
    DeformBasisType::Type{<:DeformBasis{C}}=StdDeformBasis{C};
    special_return::Type{T}=Nothing,
) where {C <: RingElem, T <: Union{Nothing, SMat}}
    @vprintln :PBWDeformations 1 "Computing Deform Basis"
    deform_basis = DeformBasisType(sp, degs)
    return all_pbwdeformations(sp, deform_basis; special_return)
end

"""
    all_pbwdeformations(sp::SmashProductLie{C}, deg::Int, DeformBasisType::Type{<:DeformBasis{C}}=StdDeformBasis{C}; special_return=Nothing) where {C <: RingElem}

The same as the other method, but only for a single degree `deg`.
"""
function all_pbwdeformations(
    sp::SmashProductLie{C},
    deg::Int,
    DeformBasisType::Type{<:DeformBasis{C}}=StdDeformBasis{C};
    special_return::Type{T}=Nothing,
) where {C <: RingElem, T <: Union{Nothing, SMat}}
    return all_pbwdeformations(sp, [deg], DeformBasisType; special_return)
end
