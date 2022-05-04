mutable struct SmashProductDeformLie{C <: RingElement}
    dimL::Int64
    dimV::Int64
    baseL::Vector{QuadraticQuoAlgebraElem{C}}
    baseV::Vector{QuadraticQuoAlgebraElem{C}}
    coeff_ring::Ring
    alg::QuadraticQuoAlgebra{C}
    # dynkin::Char
    # n::Int64
    # lambda::Vector{Int64}
    # matrixRepL::Vector{Matrix{Int64}}
    symmetric::Bool
    kappa::Matrix{QuadraticQuoAlgebraElem{C}}
end


function smash_product_deform_lie(
    sp::SmashProductLie{C},
    kappa::Matrix{QuadraticQuoAlgebraElem{C}},
) where {C <: RingElement}
    size(kappa) == (sp.dimV, sp.dimV) || throw(ArgumentError("kappa has wrong dimensions."))

    dimL = sp.dimL
    dimV = sp.dimV
    coeff_ring = sp.coeff_ring
    baseL = sp.baseL
    baseV = sp.baseV

    for i in 1:dimV, j in 1:i
        kappa[i, j] == -kappa[j, i] || throw(ArgumentError("kappa is not skew-symmetric."))
        all(x -> x <= dimL, var_ids(kappa[i, j])) ||
            throw(ArgumentError("kappa does not only take values in the hopf algebra"))
        all(x -> x <= dimL, var_ids(kappa[j, i])) ||
            throw(ArgumentError("kappa does not only take values in the hopf algebra"))
    end

    symmetric = true
    rels = Dict{Tuple{Int, Int}, QuadraticQuoAlgebraElem{C}}()
    for i in 1:dimV, j in 1:dimV
        # We have the commutator relation [v_i, v_j] = kappa[i,j]
        # which is equivalent to v_i*v_j = v_j*v_i + kappa[i,j]
        rels[(dimL + i, dimL + j)] = baseV[j] * baseV[i] + kappa[i, j]
        symmetric &= iszero(kappa[i, j])
    end

    alg, _ = quadratic_quo_algebra(sp.alg, rels)
    baseL = map(alg, baseL)
    baseV = map(alg, baseV)

    return SmashProductDeformLie{C}(dimL, dimV, baseL, baseV, coeff_ring, alg, symmetric, kappa), (baseL, baseV)
end

function smash_product_symmdeform_lie(sp::SmashProductLie{C}) where {C <: RingElement}
    kappa = fill(zero(sp.alg), sp.dimV, sp.dimV)
    return smash_product_deform_lie(sp, kappa)
end

ngens(d::SmashProductDeformLie) = d.dimL, d.dimV

function gens(d::SmashProductDeformLie{C}) where {C <: RingElement}
    return [gen(d.alg, i) for i in 1:d.dimL], [gen(d.alg, i + d.dimL) for i in 1:d.dimV]
end


function show(io::IO, deform::SmashProductDeformLie)
    local max_gens = 4 # largest number of generators to print
    if deform.symmetric
        print(io, "Symmetric ")
    end
    print(io, "Deformation of ")
    print(io, "Lie Algebra Smash Product with basis ")
    for i in 1:min(deform.dimL - 1, max_gens - 1)
        print(io, string(deform.alg.S[i]), ", ")
    end
    if deform.dimL > max_gens
        print(io, "..., ")
    end
    print(io, string(deform.alg.S[deform.dimL]) * ", ")
    for i in 1:min(deform.dimV - 1, max_gens - 1)
        print(io, string(deform.alg.S[deform.dimL+i]), ", ")
    end
    if deform.dimV > max_gens
        print(io, "..., ")
    end
    print(io, string(deform.alg.S[deform.dimL+deform.dimV]))
    print(io, " over ")
    print(IOContext(io, :compact => true), deform.coeff_ring)
end

function change_base_ring(R::Ring, d::SmashProductDeformLie{C}) where {C <: RingElement}
    alg = change_base_ring(R, d.alg)
    baseL = [gen(alg, i) for i in 1:d.dimL]
    baseV = [gen(alg, d.dimL + i) for i in 1:d.dimV]
    kappa = map(alg, d.kappa)

    return SmashProductDeformLie{elem_type(R)}(d.dimL, d.dimV, baseL, baseV, R, alg, d.symmetric, kappa)
end


function pbwdeform_eqs(deform::SmashProductDeformLie{C}; disabled::Vector{Symbol}=Symbol[]) where {C <: RingElement}
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras.
    # DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011
    dimL = deform.dimL
    dimV = deform.dimV
    kappa = deform.kappa
    x(i) = gen(deform.alg, i)
    v(i) = gen(deform.alg, dimL + i)

    ## (a) κ is H-invariant
    iter_a =
        :a in disabled ? [] :
        (
            comm(comm(h, v(i); strict=true), v(j)) # κ([h⋅v_i,v_j])
            + comm(v(i), comm(h, v(j); strict=true)) # κ([v_i,h⋅v_j])
            - comm(h, comm(v(i), v(j); strict=true)) # h⋅κ([v_i,v_j])
            for (h, (i, j)) in Iterators.product([x(i) for i in 1:dimL], Combinatorics.Combinations(dimV, 2))
        )
    # m[1][2] denotes the index of the only basis element in the monomial m

    ## (b) trivial
    iter_b = :b in disabled ? [] : []

    ## (c) 0 = κ ⊗ id - id ⊗ κ on (I ⊗ V) ∩ (V ⊗ I)
    # (I ⊗ V) ∩ (V ⊗ I) has basis v_iv_jv_k + v_jv_kv_i + v_kv_iv_j - v_kv_jv_i - v_jv_iv_k - v_iv_kv_j for i<j<k
    iter_c =
        :c in disabled ? [] :
        (
            comm(v(i), v(j); strict=true) * v(k) - v(i) * comm(v(j), v(k); strict=true) +
            comm(v(j), v(k); strict=true) * v(i) - v(j) * comm(v(k), v(i); strict=true) +
            comm(v(k), v(i); strict=true) * v(j) - v(k) * comm(v(i), v(j); strict=true) -
            comm(v(k), v(j); strict=true) * v(i) + v(k) * comm(v(j), v(i); strict=true) -
            comm(v(j), v(i); strict=true) * v(k) + v(j) * comm(v(i), v(k); strict=true) -
            comm(v(i), v(k); strict=true) * v(j) + v(i) * comm(v(k), v(j); strict=true) for
            (i, j, k) in Combinatorics.Combinations(dimV, 3)
        )

    ## (d) trivial
    iter_d = :d in disabled ? [] : []

    iter = Iterators.flatten([iter_a, iter_b, iter_c, iter_d])
    return Iterators.map(normal_form, iter)
end

function pbwdeform_neqs(deform::SmashProductDeformLie{C}) where {C <: RingElement}
    num_a = deform.dimL * binomial(deform.dimV, 2)
    num_b = 0
    num_c = binomial(deform.dimV, 3)
    num_d = 0

    return num_a + num_b + num_c + num_d
end

function ispbwdeform(d::SmashProductDeformLie{C})::Bool where {C <: RingElement}
    return all(iszero, pbwdeform_eqs(d))
end


###############################################################################
#
#   All PBW deformations
#
###############################################################################

@inline function coefficient_comparison(eq::AlgebraElem{C}) where {C <: RingElement}
    return eq.coeffs
end

@inline function linpoly_to_spvector(a::fmpq_mpoly, var_lookup::Dict{fmpq_mpoly, Int64}, nvars::Int64)
    @assert total_degree(a) <= 1

    return sparsevec(Dict(var_lookup[monomial(a, i)] => coeff(a, i) for i in 1:length(a)), nvars)
end

function reduce_and_store!(
    lgs::Vector{Union{Nothing, SparseVector{T, Int64}}},
    v::SparseVector{T, Int64},
) where {T <: Union{RingElement, Number}}
    while !iszero(v)
        nz_inds, nz_vals = findnz(v)
        v = inv(nz_vals[1]) .* v
        if lgs[nz_inds[1]] === nothing
            lgs[nz_inds[1]] = v
            return
        else
            v -= lgs[nz_inds[1]]
        end
    end
end

function reduced_row_echelon!(
    lgs::Vector{Union{Nothing, SparseVector{T, Int64}}},
) where {T <: Union{RingElement, Number}}
    for i in length(lgs):-1:1
        if lgs[i] === nothing
            continue
        end
        nz_inds, nz_vals = findnz(lgs[i])
        for (ind, j) in enumerate(nz_inds[2:end])
            if lgs[j] !== nothing
                lgs[i] -= nz_vals[ind+1] .* lgs[j]
            end
        end
    end
    return lgs
end

function lgs_to_mat(lgs::Vector{Union{Nothing, SparseVector{T, Int64}}}) where {T <: Union{RingElement, Number}}
    n = length(lgs)
    mat = spzeros(T, n, n)
    for i in 1:n
        if lgs[i] !== nothing
            mat[i, :] = lgs[i]
        end
    end
    return mat
end

function indices_of_freedom(mat::SparseArrays.SparseMatrixCSC{T, Int64}) where {T <: Union{RingElement, Number}}
    size(mat)[1] == size(mat)[2] || throw(ArgumentError("Matrix needs to be square."))
    return filter(i -> iszero(mat[i, i]), 1:size(mat)[1])
end


abstract type DeformBase{C <: RingElement} end
struct DeformStdBase{C <: RingElement} <: DeformBase{C}
    length::Int
    generator

    function DeformStdBase{C}(sp::SmashProductLie{C}, maxdeg::Int) where {C <: RingElement}
        dimL = sp.dimL
        dimV = sp.dimV
        R = coefficient_ring(sp.alg)
        generator = (
            begin
                kappa = fill(sp.alg(0), dimV, dimV)
                entry = prod(map(k -> sp.baseL[k], ind); init=sp.alg(1))
                kappa[i, j] += entry
                kappa[j, i] -= entry
                kappa
            end for i in 1:dimV for j in i+1:dimV for d in 0:maxdeg for
            ind in Combinatorics.with_replacement_combinations(1:dimL, d)
        )

        length = div(dimV * (dimV - 1), 2) * sum(binomial(dimL + k - 1, k) for k in 0:maxdeg)
        return new{C}(length, generator)
    end
end

function Base.length(base::DeformStdBase)
    return base.length
end


function pbwdeforms_all(
    sp::SmashProductLie{C},
    maxdeg::Int,
    DeformBaseType::Type{<:DeformBase{C}}=DeformStdBase{C};
    special_return::Type{T}=Nothing,
) where {C <: RingElement, T <: Union{Nothing, SparseMatrixCSC}}
    dimL = sp.dimL
    dimV = sp.dimV

    deform_base = DeformBaseType(sp, maxdeg)
    nvars = length(deform_base)

    @info "Constructing MPolyRing..."
    R, vars = PolynomialRing(sp.coeff_ring, nvars)
    var_lookup = Dict(vars[i] => i for i in 1:nvars)

    @info "Changing SmashProductLie coeffcient type..."
    new_sp = change_base_ring(R, sp)

    @info "Constructing kappa..."
    kappa = fill(new_sp.alg(0), dimV, dimV)
    for (i, b) in enumerate(deform_base.generator)
        kappa += vars[i] .* new_sp.alg.(b)
    end

    @info "Constructing deformation..."
    deform = smash_product_deform_lie(new_sp, kappa)[1]

    @info "Generating equation iterator..."
    neqs = pbwdeform_neqs(deform)
    iter = Iterators.map(
        a -> linpoly_to_spvector(a, var_lookup, nvars),
        Iterators.flatten(
            Iterators.map(
                function (x)
                    i = x[1]
                    a = x[2]
                    @debug "Equation $i/$(neqs), $(floor(Int, 100*i / neqs))%"
                    coefficient_comparison(a)
                end,
                enumerate(pbwdeform_eqs(deform)),
            ),
        ),
    )

    @info "Computing row-echelon form..."
    lgs = Vector{Union{Nothing, SparseVector{fmpq, Int64}}}(nothing, nvars)
    for v in iter
        reduce_and_store!(lgs, v)
    end

    @info "Computing reduced row-echelon form..."
    reduced_row_echelon!(lgs)

    mat = lgs_to_mat(lgs)

    if special_return === SparseMatrixCSC
        return mat, vars
    end

    @info "Computing a basis..."
    freedom_ind = indices_of_freedom(mat)
    freedom_deg = length(freedom_ind)
    kappas = Vector{Matrix{QuadraticQuoAlgebraElem{C}}}(undef, freedom_deg)
    for l in 1:freedom_deg
        kappas[l] = fill(sp.alg(0), dimV, dimV)
    end
    if freedom_deg > 0
        for (i, b) in enumerate(deform_base.generator)
            if iszero(mat[i, i])
                l = findfirst(isequal(i), freedom_ind)
                kappas[l] += b
            else
                for j in i+1:nvars
                    if !iszero(mat[i, j])
                        l = findfirst(isequal(j), freedom_ind)
                        kappas[l] += -mat[i, j] .* b
                    end
                end
            end
        end
    end
    return kappas
end
