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
        dimV, e = extract_sp_info__so_extpowers_stdmod(sp)
        extra_data = Dict{DeformationMap{C}, Set{ArcDiagram}}()
        normalize = no_normalize ? identity : normalize_default

        lens = []
        iters = []
        debug_counter = 0
        for d in degs
            diag_iter = pbw_arc_diagrams__so_extpowers_stdmod(e, d)
            len = length(diag_iter)
            iter = (
                begin
                    @debug "Basis generation deg $(d), $(debug_counter = (debug_counter % len) + 1)/$(len), $(floor(Int, 100*debug_counter / len))%"
                    basis_elem = arcdiag_to_basiselem__so_extpowers_stdmod(diag, dimV, e, d, sp.alg(0), sp.basisL)
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


function extract_sp_info__so_extpowers_stdmod(sp::SmashProductLie{C}) where {C <: RingElement}
    if isnothing(sp.info.dynkin) || isnothing(sp.info.n)
        error("Dynkin type unknown, but needed.")
    elseif !sp.info.constructive_basis
        error("Constructive basis needed.")
    elseif sp.info.dynkin == 'B'
        dimV = 2 * sp.info.n + 1
    elseif sp.info.dynkin == 'D'
        dimV = 2 * sp.info.n
    else
        error("Dynkin type '$(sp.info.dynkin)' not supported.")
    end

    if isnothing(sp.info.power_of_std_mod) || !(sp.info.power_of_std_mod < 0)
        error("Module needs to be an exterior power of the standard module.")
    end
    e = -sp.info.power_of_std_mod

    return dimV, e
end

function pbw_arc_diagrams__so_extpowers_stdmod(e::Int, d::Int)
    indep_sets = [1:e, e+1:2*e, [[2 * e + 2 * i - 1, 2 * e + 2 * i] for i in 1:d]...]
    return all_arc_diagrams(2 * e, 2 * d; indep_sets)
end


function arcdiag_to_basiselem__so_extpowers_stdmod(
    diag::ArcDiagram,
    dimV::Int,
    e::Int,
    d::Int,
    zero::QuadraticQuoAlgebraElem{C},
    basisL::Vector{QuadraticQuoAlgebraElem{C}},
) where {C <: RingElement}
    iso_wedge2V_g = Dict{Vector{Int}, Int}()
    for (i, bs) in enumerate(Combinatorics.combinations(1:dimV, 2))
        iso_wedge2V_g[bs] = i
    end
    index = Dict{Vector{Int}, Int}()
    for (i, is) in enumerate(Combinatorics.combinations(1:dimV, e))
        index[is] = i
    end

    kappa = fill(zero, binomial(dimV, e), binomial(dimV, e))
    for is in Combinatorics.combinations(1:dimV, e), js in Combinatorics.combinations(1:dimV, e)
        i = index[is]
        j = index[js]
        if i >= j
            continue
        end
        for is in Combinatorics.permutations(is), js in Combinatorics.permutations(js), swap in [false, true]
            sgn_upper_labels = levicivita(sortperm(is)) * levicivita(sortperm(js))
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
            free_index = Dict{Int, Int}()
            for (k, f) in enumerate(frees)
                free_index[f] = k
            end
            entry = zero
            for labeling in (isempty(frees) ? [Int[]] : AbstractAlgebra.ProductIterator(1:dimV, length(frees)))
                lower_labeled = [
                    labeled_diag[k] != 0 ? labeled_diag[k] : labeling[free_index[min(k, diag.adjacency[k])]] for
                    k in 2*e+1:2*e+2*d
                ]
                zeroelem = false
                sign_lower_labels = 1
                basiselem = Int[]
                for k in 1:2:length(lower_labeled)
                    if lower_labeled[k] == lower_labeled[k+1]
                        zeroelem = true
                        break
                    elseif lower_labeled[k] > lower_labeled[k+1]
                        sign_lower_labels *= -1
                        append!(basiselem, iso_wedge2V_g[[lower_labeled[k+1], lower_labeled[k]]])
                    else
                        append!(basiselem, iso_wedge2V_g[[lower_labeled[k], lower_labeled[k+1]]])
                    end
                end
                if zeroelem
                    continue
                end
                symm_basiselem =
                    1 // factorial(length(basiselem)) *
                    sum(prod(basisL[ind]) for ind in Combinatorics.permutations(basiselem))
                entry += sign_lower_labels * normal_form(symm_basiselem)
            end

            entry *= sgn_upper_labels

            kappa[i, j] += entry
            kappa[j, i] -= entry
        end
    end
    return kappa
end
