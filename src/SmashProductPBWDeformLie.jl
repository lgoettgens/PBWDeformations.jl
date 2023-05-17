"""
    pbwdeform_eqs(d::SmashProductDeformLie{C}; disabled::Vector{Symbol}=Symbol[]) where {C <: RingElement}

Returns the equations for `d` being a PBW deformation of a smash product as
in Theorem 3.1 of [WW14](@cite).
Subsets of the equations can be disabled by passing the corresponding symbols as
keyword arguments, e.g. `disabled = [:c, :d]`.
"""
function pbwdeform_eqs(d::SmashProductDeformLie{C}; disabled::Vector{Symbol}=Symbol[]) where {C <: RingElement}
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras.
    # DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011
    dimL = dim(d.sp.L)
    dimV = dim(d.sp.V)

    x(i) = gen(d, i, :L)
    v(i) = gen(d, i, :V)

    nfcomm(a, b) = normal_form(comm(a, b), d.rels)

    ## (a) κ is H-invariant
    iter_a =
        :a in disabled ? [] :
        (
            comm(nfcomm(h, v(i)), v(j)) # κ([h⋅v_i,v_j])
            + comm(v(i), nfcomm(h, v(j))) # κ([v_i,h⋅v_j])
            - comm(h, nfcomm(v(i), v(j))) # h⋅κ([v_i,v_j])
            for (h, (i, j)) in Iterators.product([x(i) for i in 1:dimL], combinations(dimV, 2))
        )

    ## (b) trivial
    iter_b = :b in disabled ? [] : []

    ## (c) 0 = κ ⊗ id - id ⊗ κ on (I ⊗ V) ∩ (V ⊗ I)
    # (I ⊗ V) ∩ (V ⊗ I) has basis v_iv_jv_k + v_jv_kv_i + v_kv_iv_j - v_kv_jv_i - v_jv_iv_k - v_iv_kv_j for i<j<k
    iter_c =
        :c in disabled ? [] :
        (
            nfcomm(v(i), v(j)) * v(k) - v(i) * nfcomm(v(j), v(k)) + nfcomm(v(j), v(k)) * v(i) -
            v(j) * nfcomm(v(k), v(i)) + nfcomm(v(k), v(i)) * v(j) - v(k) * nfcomm(v(i), v(j)) -
            nfcomm(v(k), v(j)) * v(i) + v(k) * nfcomm(v(j), v(i)) - nfcomm(v(j), v(i)) * v(k) +
            v(j) * nfcomm(v(i), v(k)) - nfcomm(v(i), v(k)) * v(j) + v(i) * nfcomm(v(k), v(j)) for
            (i, j, k) in combinations(dimV, 3)
        )

    ## (d) trivial
    iter_d = :d in disabled ? [] : []

    iter = Iterators.flatten([iter_a, iter_b, iter_c, iter_d])
    return Iterators.map(x -> normal_form(x, d.rels), iter)
end

function pbwdeform_neqs(d::SmashProductDeformLie{C}) where {C <: RingElement}
    dimL = dim(d.sp.L)
    dimV = dim(d.sp.V)

    num_a = dimL * binomial(dimV, 2)
    num_b = 0
    num_c = binomial(dimV, 3)
    num_d = 0

    return num_a + num_b + num_c + num_d
end

"""
    is_pbwdeformation(d::SmashProductDeformLie{C}) where {C <: RingElement}

Check if `d` is a Poincare-Birkhoff-Witt deformation of a smash product.
Uses Theorem 3.1 of [WW14](@cite).
"""
function is_pbwdeformation(d::SmashProductDeformLie{C}) where {C <: RingElement}
    return all(iszero, pbwdeform_eqs(d))
end


###############################################################################
#
#   All PBW deformations
#
###############################################################################

@inline function coefficient_comparison(eq::FreeAssAlgElem{C}) where {C <: RingElement}
    return collect(coefficients(eq))
end

@inline function linpoly_to_spvector(a::QQMPolyRingElem, var_lookup::Dict{QQMPolyRingElem, Int})
    @assert total_degree(a) <= 1

    return sparse_row(QQ, [(var_lookup[monomial(a, i)], coeff(a, i)) for i in 1:length(a)])
end

function reduce_and_store!(lgs::Vector{Union{Nothing, SRow{T}}}, v::SRow{T}) where {T <: RingElement}
    while !is_zero(v)
        nz_ind, nz_val = first(v)

        v = divexact(v, nz_val)

        if isnothing(lgs[nz_ind])
            lgs[nz_ind] = v
            return
        else
            v -= lgs[nz_ind]
        end
    end
end

function reduced_row_echelon!(lgs::Vector{Union{Nothing, SRow{T}}}) where {T <: RingElement}
    for i in length(lgs):-1:1
        if isnothing(lgs[i])
            continue
        end
        for (nz_ind, nz_val) in Iterators.drop(lgs[i], 1)
            if !isnothing(lgs[nz_ind])
                lgs[i] -= nz_val * lgs[nz_ind]
            end
        end
    end
    return lgs
end

function lgs_to_mat(lgs::Vector{Union{Nothing, SRow{T}}}, R::Ring) where {T <: RingElement}
    n = length(lgs)
    mat = sparse_matrix(R, n, n)
    for i in 1:n
        if !isnothing(lgs[i])
            setindex!(mat, lgs[i], i)
        end
    end
    return mat
end

function indices_of_freedom(mat::SMat{T}) where {T <: RingElement}
    return map(j -> (Oscar.Hecke.find_row_starting_with(mat, j), j), 1:size(mat)[2]) |>
           filter(ij -> is_zero(mat[ij[1], ij[2]])) |>
           x -> map(ij -> ij[2], x)
end


"""
    all_pbwdeformations(sp::SmashProductLie{C}, deform_basis::DeformBasis{C}; special_return=Nothing) where {C <: RingElement}

Computes a basis of all Poincare-Birkhoff-Witt deformations of `sp`.
`deform_basis` specifies the basis to use for the space of deformation maps.
If `special_return` is `SMat`, the function returns intermediate results.

Uses [`pbwdeform_eqs`](@ref) and thus Theorem 3.1 of [WW14](@cite).
"""
function all_pbwdeformations(
    sp::SmashProductLie{C},
    deform_basis::DeformBasis{C};
    special_return::Type{T}=Nothing,
) where {C <: RingElement, T <: Union{Nothing, SMat}}
    @req sp.coeff_ring == QQ "Only implemented for QQ coefficients."

    dimL = dim(sp.L)
    dimV = dim(sp.V)

    nvars = length(deform_basis)

    @info "Constructing MPolyRing..."
    R, vars = polynomial_ring(sp.coeff_ring, max(nvars, 1))
    var_lookup = Dict(vars[i] => i for i in 1:nvars)

    @info "Changing SmashProductLie coeffcient type..."
    new_sp = change_base_ring(R, sp)

    @info "Constructing kappa..."
    kappa = fill(new_sp.alg(0), dimV, dimV)
    for (i, b) in enumerate(deform_basis)
        kappa += vars[i] .* map(e -> change_base_ring(R, e, parent=new_sp.alg), b)
    end

    @info "Constructing deformation..."
    d = deform(new_sp, kappa)

    @info "Generating equation iterator..."
    neqs = pbwdeform_neqs(d)
    iter = Iterators.map(
        a -> linpoly_to_spvector(a, var_lookup),
        Iterators.flatten(
            Iterators.map(function (x)
                    i = x[1]
                    a = x[2]
                    @debug "Equation $i/$(neqs), $(floor(Int, 100*i / neqs))%"
                    coefficient_comparison(a)
                end, enumerate(pbwdeform_eqs(d))),
        ),
    )

    @info "Computing matrix..."
    lgs = sparse_matrix(sp.coeff_ring, 0, nvars)
    for v in iter
        push!(lgs, v)
    end

    @info "Computing reduced row-echelon form..."
    rref!(lgs; truncate=true)

    if special_return <: SMat
        return lgs, vars
    end

    @info "Computing a basis..."
    freedom_ind = indices_of_freedom(lgs)
    freedom_deg = length(freedom_ind)
    kappas = Vector{DeformationMap{C}}(undef, freedom_deg)
    for l in 1:freedom_deg
        kappas[l] = fill(sp.alg(0), dimV, dimV)
    end
    if freedom_deg > 0
        for (j, b) in enumerate(deform_basis)
            l = findfirst(==(j), freedom_ind)
            if !isnothing(l)
                kappas[l] += b
            else
                i = Oscar.Hecke.find_row_starting_with(lgs, j)
                for k in j+1:nvars
                    if !iszero(lgs[i, k])
                        l = findfirst(==(k), freedom_ind)
                        kappas[l] += -lgs[i, k] .* b
                    end
                end
            end
        end
    end
    return kappas
end

"""
    all_pbwdeformations(sp::SmashProductLie{C}, degs::AbstractVector{Int}, DeformBasisType::Type{<:DeformBasis{C}}=StdDeformBasis{C}; special_return=Nothing) where {C <: RingElement}

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
) where {C <: RingElement, T <: Union{Nothing, SMat}}
    @info "Computing Deform Basis"
    deform_basis = DeformBasisType(sp, degs)
    return all_pbwdeformations(sp, deform_basis; special_return)
end

"""
    all_pbwdeformations(sp::SmashProductLie{C}, deg::Int, DeformBasisType::Type{<:DeformBasis{C}}=StdDeformBasis{C}; special_return=Nothing) where {C <: RingElement}

The same as the other method, but only for a single degree `deg`.
"""
function all_pbwdeformations(
    sp::SmashProductLie{C},
    deg::Int,
    DeformBasisType::Type{<:DeformBasis{C}}=StdDeformBasis{C};
    special_return::Type{T}=Nothing,
) where {C <: RingElement, T <: Union{Nothing, SMat}}
    return all_pbwdeformations(sp, [deg], DeformBasisType; special_return)
end
