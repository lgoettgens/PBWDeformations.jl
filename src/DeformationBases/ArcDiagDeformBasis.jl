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
        get_attribute(sp.L, :type, nothing) == :special_orthogonal || error("Only works for so_n.")
        sp.V isa Union{LieAlgebraExteriorPowerModule{C}, LieAlgebraSymmetricPowerModule{C}} &&
            sp.V.inner_mod isa LieAlgebraStdModule{C} || error("Only works for exterior powers of the standard module.")

        extra_data = Dict{DeformationMap{C}, Set{ArcDiagram}}()
        normalize = no_normalize ? identity : normalize_default

        lens = []
        iters = []
        debug_counter = 0
        for d in degs
            diag_iter = pbw_arc_diagrams__so(sp.V, d)
            len = length(diag_iter)
            iter = (
                begin
                    @debug "Basis generation deg $(d), $(debug_counter = (debug_counter % len) + 1)/$(len), $(floor(Int, 100*debug_counter / len))%"
                    basis_elem = arcdiag_to_deformationmap__so(diag, sp)
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


function pbw_arc_diagrams__so(V::LieAlgebraModule{C}, d::Int) where {C <: RingElement}
    e = arc_diagram_num_points__so(V)
    upper_indep_sets = Vector{Int}[is .+ a * e for a in [0, 1] for is in arc_diagram_indep_sets__so(V)]
    lower_indep_sets = Vector{Int}[[[2i - 1, 2i] for i in 1:d]...]
    indep_sets = Vector{Int}[upper_indep_sets..., [is .+ 2e for is in lower_indep_sets]...]
    return all_arc_diagrams(2e, 2d; indep_sets)
end

function arc_diagram_num_points__so(_::LieAlgebraStdModule{C}) where {C <: RingElement}
    return 1
end

function arc_diagram_num_points__so(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement}
    return arc_diagram_num_points__so(V.inner_mod) * V.power
end

function arc_diagram_num_points__so(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement}
    return arc_diagram_num_points__so(V.inner_mod) * V.power
end

function arc_diagram_num_points__so(V::LieAlgebraTensorPowerModule{C}) where {C <: RingElement}
    return arc_diagram_num_points__so(V.inner_mod) * V.power
end


function arc_diagram_indep_sets__so(_::LieAlgebraStdModule{C}) where {C <: RingElement}
    return Vector{Int}[]
end

function arc_diagram_indep_sets__so(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement}
    if V.inner_mod isa LieAlgebraStdModule{C}
        return [1:V.power]
    else
        is = arc_diagram_indep_sets__so(V.inner_mod)
        return [map(i -> i + k * arc_diagram_num_points__so(V.inner_mod), is) for k in 0:V.power-1, is in is]
    end
end

function arc_diagram_indep_sets__so(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement}
    if V.inner_mod isa LieAlgebraStdModule{C}
        return Vector{Int}[]
    else
        is = arc_diagram_indep_sets__so(V.inner_mod)
        return [map(i -> i + k * arc_diagram_num_points__so(V.inner_mod), is) for k in 0:V.power-1, is in is]
    end
end

function arc_diagram_indep_sets__so(V::LieAlgebraTensorPowerModule{C}) where {C <: RingElement}
    if V.inner_mod isa LieAlgebraStdModule{C}
        return Vector{Int}[]
    else
        is = arc_diagram_indep_sets__so(V.inner_mod)
        return [map(i -> i + k * arc_diagram_num_points__so(V.inner_mod), is) for k in 0:V.power-1, is in is]
    end
end


function arc_diagram_label_iterator__so(_::LieAlgebraStdModule, base_labels::AbstractVector{Int})
    return [[l] for l in base_labels]
end

function arc_diagram_label_iterator__so(V::LieAlgebraExteriorPowerModule, base_labels::AbstractVector{Int})
    return Combinatorics.combinations(collect(arc_diagram_label_iterator__so(V.inner_mod, base_labels)), V.power) .|>
           Iterators.flatten .|>
           collect

end

function arc_diagram_label_iterator__so(V::LieAlgebraSymmetricPowerModule, base_labels::AbstractVector{Int})
    return Combinatorics.with_replacement_combinations(
               collect(arc_diagram_label_iterator__so(V.inner_mod, base_labels)),
               V.power,
           ) .|>
           Iterators.flatten .|>
           collect
end

function arc_diagram_label_iterator__so(V::LieAlgebraTensorPowerModule, base_labels::AbstractVector{Int})
    return ProductIterator(arc_diagram_label_iterator__so(V.inner_mod, base_labels), V.power) .|>
           reverse .|>
           Iterators.flatten .|>
           collect
end


function arc_diagram_label_permutations__so(_::LieAlgebraStdModule, label::AbstractVector{Int})
    length(label) == 1 || error("Number of labels mistmatch.")
    return [(label, 1)]
end

function arc_diagram_label_permutations__so(V::LieAlgebraExteriorPowerModule, label::AbstractVector)
    m = arc_diagram_num_points__so(V.inner_mod)
    length(label) == m * V.power || error("Number of labels mistmatch.")
    return [
        begin
            inner_label = flatten(first.(inner_iter))
            inner_sign = prod(last.(inner_iter))
            (inner_label, inner_sign * levicivita(outer_perm))
        end for outer_perm in Combinatorics.permutations(1:V.power) for inner_iter in ProductIterator([
            arc_diagram_label_permutations__so(V.inner_mod, label[(outer_perm[i]-1)*m+1:outer_perm[i]*m]) for
            i in 1:V.power
        ])
    ]
end

function arc_diagram_label_permutations__so(V::LieAlgebraSymmetricPowerModule, label::AbstractVector)
    m = arc_diagram_num_points__so(V.inner_mod)
    length(label) == m * V.power || error("Number of labels mistmatch.")
    return [
        begin
            inner_label = flatten(first.(inner_iter))
            inner_sign = prod(last.(inner_iter))
            (inner_label, inner_sign)
        end for outer_perm in Combinatorics.permutations(1:V.power) for inner_iter in ProductIterator([
            arc_diagram_label_permutations__so(V.inner_mod, label[(outer_perm[i]-1)*m+1:outer_perm[i]*m]) for
            i in 1:V.power
        ])
    ]
end

function arc_diagram_label_permutations__so(V::LieAlgebraTensorPowerModule, label::AbstractVector)
    m = arc_diagram_num_points__so(V.inner_mod)
    length(label) == m * V.power || error("Number of labels mistmatch.")
    return [
        begin
            inner_label = flatten(first.(inner_iter))
            inner_sign = prod(last.(inner_iter))
            (inner_label, inner_sign)
        end for inner_iter in
        ProductIterator([arc_diagram_label_permutations__so(V.inner_mod, label[(i-1)*m+1:i*m]) for i in 1:V.power])
    ]
end


function arcdiag_to_deformationmap__so(diag::ArcDiagram, sp::SmashProductLie{C}) where {C <: RingElement}
    d = div(diag.nLower, 2)
    dim_stdmod_V = sp.L.n

    e = arc_diagram_num_points__so(sp.V)

    iso_wedge2V_g = Dict{Vector{Int}, Int}()
    for (i, bs) in enumerate(Combinatorics.combinations(1:sp.L.n, 2))
        iso_wedge2V_g[bs] = i
    end

    index = Dict{Vector{Int}, Int}()
    for (i, is) in enumerate(arc_diagram_label_iterator__so(sp.V, 1:dim_stdmod_V))
        index[is] = i
    end

    kappa = fill(zero(sp.alg), dim(sp.V), dim(sp.V))
    for is in arc_diagram_label_iterator__so(sp.V, 1:dim_stdmod_V),
        js in arc_diagram_label_iterator__so(sp.V, 1:dim_stdmod_V)

        i = index[is]
        j = index[js]
        if i >= j
            continue
        end
        for (is, sgn_left) in arc_diagram_label_permutations__so(sp.V, is),
            (js, sgn_right) in arc_diagram_label_permutations__so(sp.V, js),
            swap in [false, true]

            sgn_upper_labels = sgn_left * sgn_right
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
            entry = zero(sp.alg)

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
                        symm_basiselem = sp.alg(
                            fill(C(1 // factorial(length(basiselem))), factorial(length(basiselem))),
                            [ind for ind in Combinatorics.permutations(basiselem)],
                        )
                        entry += sign_lower_labels * normal_form(symm_basiselem, sp.rels)
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
