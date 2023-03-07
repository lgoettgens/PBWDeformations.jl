"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by an arc diagram of a suitable size,
which gets symmetrized and specialised to the given smash product.
This process is due to [FM22](@cite).
"""
struct ArcDiagDeformBasis{C <: RingElement} <: DeformBasis{C}
    len::Int
    iter
    extra_data::Dict{DeformationMap{C}, Set{ArcDiagram}}
    normalize

    function ArcDiagDeformBasis{C}(
        sp::SmashProductLie{C},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElement}
        n, e, typeof_power = extract_sp_info__so_powers_stdmod(sp)
        extra_data = Dict{DeformationMap{C}, Set{ArcDiagram}}()
        normalize = no_normalize ? identity : normalize_default

        lens = []
        iters = []
        debug_counter = 0
        for d in degs
            diag_iter = pbw_arc_diagrams__so_powers_stdmod(typeof_power, e, d)
            len = length(diag_iter)
            iter = (
                begin
                    @debug "Basis generation deg $(d), $(debug_counter = (debug_counter % len) + 1)/$(len), $(floor(Int, 100*debug_counter / len))%"
                    basis_elem =
                        arcdiag_to_basiselem__so_powers_stdmod(diag, n, typeof_power, e, d, sp.alg(0), sp.rels)
                    if !no_normalize
                        basis_elem = normalize(basis_elem)
                    end
                    if haskey(extra_data, basis_elem)
                        push!(extra_data[basis_elem], diag)
                    else
                        extra_data[basis_elem] = Set([diag])
                    end
                    basis_elem
                end for
                diag in diag_iter if is_crossing_free(diag, part=:upper) && is_crossing_free(diag, part=:lower)
            )
            push!(lens, len)
            push!(iters, iter)
        end
        len = sum(lens)
        iter = Iterators.flatten(iters)
        if !no_normalize
            iter = unique(Iterators.filter(b -> !iszero(b), iter))
            len = length(iter)
        end
        return new{C}(len, iter, extra_data, normalize)
    end
end

function Base.iterate(i::ArcDiagDeformBasis)
    return iterate(i.iter)
end

function Base.iterate(i::ArcDiagDeformBasis, s)
    return iterate(i.iter, s)
end

Base.length(basis::ArcDiagDeformBasis) = basis.len


function extract_sp_info__so_powers_stdmod(sp::SmashProductLie{C}) where {C <: RingElement}
    has_attribute(sp, :base_liealgebra) || error("Metadata not found, but needed.")
    has_attribute(sp, :base_liealgebra_module) || error("Metadata not found, but needed.")

    L = get_attribute(sp, :base_liealgebra)
    V = get_attribute(sp, :base_liealgebra_module)

    get_attribute(L, :type) == :special_orthogonal || error("Only implemented for so_n.")
    n = L.n

    if V isa LieAlgebraExteriorPowerModule{C} && V.inner_mod isa LieAlgebraStdModule{C}
        e = V.power
        typeof_power = :exterior
    elseif V isa LieAlgebraSymmetricPowerModule{C} && V.inner_mod isa LieAlgebraStdModule{C}
        e = V.power
        typeof_power = :symmetric
    else
        error("Module needs to be an exterior or symmetric power of the standard module.")
    end

    return n, e, typeof_power
end

function pbw_arc_diagrams__so_powers_stdmod(typeof_power::Symbol, e::Int, d::Int)
    if typeof_power == :exterior
        indep_sets = [1:e, e+1:2*e, [[2 * e + 2 * i - 1, 2 * e + 2 * i] for i in 1:d]...]
        return all_arc_diagrams(2 * e, 2 * d; indep_sets)
    elseif typeof_power == :symmetric
        indep_sets = Vector{Vector{Int}}([[[2 * e + 2 * i - 1, 2 * e + 2 * i] for i in 1:d]...])
        return all_arc_diagrams(2 * e, 2 * d; indep_sets)
    else
        error("Unknown type of power.")
    end
end


function arcdiag_to_basiselem__so_powers_stdmod(
    diag::ArcDiagram,
    dim_stdmod_V::Int,
    typeof_power::Symbol,
    e::Int,
    d::Int,
    zero::FreeAssAlgElem{C},
    rels::QuadraticRelations{C},
) where {C <: RingElement}
    iso_wedge2V_g = Dict{Vector{Int}, Int}()
    for (i, bs) in enumerate(Combinatorics.combinations(1:dim_stdmod_V, 2))
        iso_wedge2V_g[bs] = i
    end

    if typeof_power == :exterior
        upper_label_iterator = Combinatorics.combinations(1:dim_stdmod_V, e)
    elseif typeof_power == :symmetric
        upper_label_iterator = Combinatorics.with_replacement_combinations(1:dim_stdmod_V, e)
    else
        error("Unknown type of power.")
    end
    index = Dict{Vector{Int}, Int}()
    for (i, is) in enumerate(upper_label_iterator)
        index[is] = i
    end

    if typeof_power == :exterior
        kappadim = binomial(dim_stdmod_V, e)
    elseif typeof_power == :symmetric
        kappadim = binomial(dim_stdmod_V + e - 1, e)
    else
        error("Unknown type of power.")
    end
    kappa = fill(zero, kappadim, kappadim)
    for is in upper_label_iterator, js in upper_label_iterator
        i = index[is]
        j = index[js]
        if i >= j
            continue
        end
        for is in Combinatorics.permutations(is), js in Combinatorics.permutations(js), swap in [false, true]
            if typeof_power == :exterior
                sgn_upper_labels = levicivita(sortperm(is)) * levicivita(sortperm(js))
            elseif typeof_power == :symmetric
                sgn_upper_labels = 1
            else
                error("Unknown type of power.")
            end
            if swap
                is, js = js, is
                sgn_upper_labels *= -1
            end

            zeroprod = false
            labeled_diag = [is..., js..., [0 for _ in 1:2d]...]
            frees = Int[]
            for k in 1:length(labeled_diag)
                if labeled_diag[k] != 0
                    if labeled_diag[diag.adjacency[k]] == 0
                        labeled_diag[diag.adjacency[k]] = labeled_diag[k]
                    else
                        if labeled_diag[k] != labeled_diag[diag.adjacency[k]]
                            zeroprod = true
                            break
                        end
                    end
                else
                    if labeled_diag[diag.adjacency[k]] == 0
                        append!(frees, min(k, diag.adjacency[k]))
                    end
                end
            end
            if zeroprod
                continue
            end
            unique!(sort!(frees))
            entry = zero

            # iterate over lower point labelings
            nextindex = 1
            while true
                if nextindex > length(frees)
                    # begin inner
                    zeroelem = false
                    sign_lower_labels = 1
                    basiselem = Int[]
                    for k in 2*e+1:2:length(labeled_diag)
                        if labeled_diag[k] == labeled_diag[k+1]
                            zeroelem = true
                            break
                        elseif labeled_diag[k] > labeled_diag[k+1]
                            sign_lower_labels *= -1
                            append!(basiselem, iso_wedge2V_g[[labeled_diag[k+1], labeled_diag[k]]])
                        else
                            append!(basiselem, iso_wedge2V_g[[labeled_diag[k], labeled_diag[k+1]]])
                        end
                    end
                    if !zeroelem
                        symm_basiselem = parent(zero)(
                            fill(C(1 // factorial(length(basiselem))), factorial(length(basiselem))),
                            [ind for ind in Combinatorics.permutations(basiselem)],
                        )
                        entry += sign_lower_labels * normal_form(symm_basiselem, rels)
                    end
                    # end inner

                    nextindex -= 1
                end

                while nextindex >= 1 && labeled_diag[frees[nextindex]] == dim_stdmod_V
                    labeled_diag[frees[nextindex]] = 0
                    labeled_diag[diag.adjacency[frees[nextindex]]] = 0
                    nextindex -= 1
                end
                if nextindex == 0
                    break
                end
                labeled_diag[frees[nextindex]] += 1
                labeled_diag[diag.adjacency[frees[nextindex]]] += 1
                if ispairgood(labeled_diag, frees[nextindex]) &&
                   ispairgood(labeled_diag, diag.adjacency[frees[nextindex]])
                    nextindex += 1
                end
            end

            entry *= sgn_upper_labels

            kappa[i, j] += entry
            kappa[j, i] -= entry
        end
    end
    return kappa
end

function ispairgood(labeled_diag::Vector{Int}, k::Int)
    left_k = k % 2 == 1 ? k : k - 1
    return labeled_diag[left_k] != labeled_diag[left_k+1]
end
