"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by an arc diagram of a suitable size,
which gets symmetrized and specialised to the given smash product.
This process is due to [FM22](@cite).
"""
struct ArcDiagDeformBasis{C <: RingElem} <: DeformBasis{C}
    len::Int
    iter
    extra_data::Dict{DeformationMap{C}, Set{ArcDiagram}}
    normalize

    function ArcDiagDeformBasis{C}(
        sp::SmashProductLie{C},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElem}
        @req get_attribute(base_lie_algebra(sp), :type, nothing) == :special_orthogonal "Only works for so_n."
        V = base_module(sp)

        V_nice, h = isomorphic_module_with_simple_structure(V)
        if !is_direct_sum(V_nice)
            temp = direct_sum(V_nice)
            h = compose(h, hom(V_nice, temp, identity_matrix(coefficient_ring(temp), dim(temp))))
            V_nice = temp
        end

        extra_data = Dict{DeformationMap{C}, Set{ArcDiagram}}()
        normalize = no_normalize ? identity : normalize_default

        lens = []
        iters = []
        debug_counter = 0
        for d in degs
            for (i_l, V_nice_summand_i_l) in enumerate(base_modules(V_nice)),
                (i_r, V_nice_summand_i_r) in enumerate(base_modules(V_nice))

                if i_l > i_r
                    continue
                end

                proj_to_summand_l = compose(h, canonical_projection(V_nice, i_l))
                proj_to_summand_r = compose(h, canonical_projection(V_nice, i_r))

                W = if i_l == i_r
                    exterior_power(V_nice_summand_i_l, 2)
                else
                    tensor_product(V_nice_summand_i_l, V_nice_summand_i_r)
                end

                diag_iter = pbw_arc_diagrams__so(W, d)
                len = length(diag_iter)
                iter = (
                    begin
                        @vprintln :PBWDeformations 2 "Basis generation deg $(lpad(d, maximum(ndigits, degs))), $(lpad(floor(Int, 100*(debug_counter = (debug_counter % len) + 1) / len), 3))%, $(lpad(debug_counter, ndigits(len)))/$(len)"
                        _basis_elem = arcdiag_to_deformationmap__so(diag, sp, W)
                        basis_elem = matrix(proj_to_summand_l) * _basis_elem * transpose(matrix(proj_to_summand_r))
                        if i_l != i_r
                            basis_elem -= transpose(basis_elem)
                        end
                        @assert is_skew_symmetric(basis_elem)

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
                #push!(iters, iter)
                push!(iters, collect(iter))
            end
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


function pbw_arc_diagrams__so(V::LieAlgebraModule, d::Int)
    n_upper_vertices = arc_diagram_num_points__so(V)
    n_lower_vertices = 2d
    upper_indep_sets = arc_diagram_indep_sets__so(V)
    lower_indep_sets = Vector{Int}[[[2i - 1, 2i] for i in 1:div(n_lower_vertices, 2)]...]
    indep_sets = Vector{Int}[[(-1) .* is for is in upper_indep_sets]; [is for is in lower_indep_sets]]
    return all_arc_diagrams(Undirected, n_upper_vertices, n_lower_vertices; indep_sets)
end

function arc_diagram_num_points__so(V::LieAlgebraModule)
    if is_standard_module(V)
        return 1
    elseif is_tensor_product(V)
        return sum(arc_diagram_num_points__so(W) for W in base_modules(V))
    elseif is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V)
        return arc_diagram_num_points__so(base_module(V)) * get_attribute(V, :power)
    else
        error("Not implemented.")
    end
end

function arc_diagram_indep_sets__so(V::LieAlgebraModule)
    if is_standard_module(V)
        return Vector{Int}[]
    elseif is_tensor_product(V)
        return Vector{Int}[] # TODO: proper implementation
    elseif is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V)
        inner_mod = base_module(V)
        power = get_attribute(V, :power)
        if is_standard_module(inner_mod)
            if is_exterior_power(V)
                return [1:power]
            else
                return Vector{Int}[]
            end
        else
            iss = arc_diagram_indep_sets__so(inner_mod)
            return [map(i -> i + k * arc_diagram_num_points__so(inner_mod), is) for k in 0:power-1 for is in iss]
        end
    else
        error("Not implemented.")
    end
end

function arc_diagram_label_iterator__so(V::LieAlgebraModule, base_labels::AbstractVector{Int})
    if is_standard_module(V)
        return [[l] for l in base_labels]
    elseif is_tensor_product(V)
        inner_mods = base_modules(V)
        return ProductIterator([
                   arc_diagram_label_iterator__so(inner_mod, base_labels) for inner_mod in reverse(inner_mods)
               ]) .|>
               reverse .|>
               Iterators.flatten .|>
               collect
    elseif is_exterior_power(V)
        inner_mod = base_module(V)
        power = get_attribute(V, :power)
        return combinations(collect(arc_diagram_label_iterator__so(inner_mod, base_labels)), power) .|>
               Iterators.flatten .|>
               collect
    elseif is_symmetric_power(V)
        inner_mod = base_module(V)
        power = get_attribute(V, :power)
        return multicombinations(collect(arc_diagram_label_iterator__so(inner_mod, base_labels)), power) .|>
               Iterators.flatten .|>
               collect
    elseif is_tensor_power(V)
        inner_mod = base_module(V)
        power = get_attribute(V, :power)
        return ProductIterator(arc_diagram_label_iterator__so(inner_mod, base_labels), power) .|>
               reverse .|>
               Iterators.flatten .|>
               collect
    else
        error("Not implemented.")
    end
end

function basis_index_mapping(V::LieAlgebraModule)
    if is_tensor_product(V)
        inner_mods = base_modules(V)
        return ProductIterator([1:dim(inner_mod) for inner_mod in reverse(inner_mods)]) .|> reverse .|> collect
    elseif is_exterior_power(V)
        inner_mod = base_module(V)
        power = get_attribute(V, :power)
        return combinations(dim(inner_mod), power) .|> collect
    elseif is_symmetric_power(V)
        inner_mod = base_module(V)
        power = get_attribute(V, :power)
        return multicombinations(dim(inner_mod), power) .|> collect
    elseif is_tensor_power(V)
        inner_mod = base_module(V)
        power = get_attribute(V, :power)
        return ProductIterator(1:dim(inner_mod), power) .|> reverse .|> collect
    else
        error("Not implemented.")
    end
end

function arc_diagram_label_permutations__so(V::LieAlgebraModule, label::AbstractVector{Int})
    if is_standard_module(V)
        @req length(label) == 1 "Number of labels mismatch."
        return [(label, 1)]
    elseif is_tensor_product(V)
        inner_mods = base_modules(V)
        @req length(label) == sum(arc_diagram_num_points__so.(inner_mods)) "Number of labels mismatch."
        return [
            begin
                inner_label = vcat(first.(inner_iter)...)
                inner_sign = prod(last.(inner_iter))
                (inner_label, inner_sign)
            end for inner_iter in ProductIterator([
                arc_diagram_label_permutations__so(
                    inner_mod,
                    label[sum(arc_diagram_num_points__so, inner_mods[1:i-1]; init=0)+1:sum(
                        arc_diagram_num_points__so,
                        inner_mods[1:i];
                        init=0,
                    )],
                ) for (i, inner_mod) in enumerate(inner_mods)
            ])
        ]
    elseif is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V)
        inner_mod = base_module(V)
        power = get_attribute(V, :power)
        m = arc_diagram_num_points__so(inner_mod)
        @req length(label) == m * power "Number of labels mismatch."
        if is_exterior_power(V)
            return [
                begin
                    inner_label = vcat(first.(inner_iter)...)
                    inner_sign = prod(last.(inner_iter))
                    (inner_label, inner_sign * outer_sign)
                end for (outer_perm, outer_sign) in permutations_with_sign(1:power) for inner_iter in ProductIterator([
                    arc_diagram_label_permutations__so(inner_mod, label[(outer_perm[i]-1)*m+1:outer_perm[i]*m]) for
                    i in 1:power
                ])
            ]
        elseif is_symmetric_power(V)
            return [
                begin
                    inner_label = vcat(first.(inner_iter)...)
                    inner_sign = prod(last.(inner_iter))
                    (inner_label, inner_sign)
                end for outer_perm in permutations(1:power) for inner_iter in ProductIterator([
                    arc_diagram_label_permutations__so(inner_mod, label[(outer_perm[i]-1)*m+1:outer_perm[i]*m]) for
                    i in 1:power
                ])
            ]
        elseif is_tensor_power(V)
            return [
                begin
                    inner_label = vcat(first.(inner_iter)...)
                    inner_sign = prod(last.(inner_iter))
                    (inner_label, inner_sign)
                end for inner_iter in
                ProductIterator([arc_diagram_label_permutations__so(inner_mod, label[(i-1)*m+1:i*m]) for i in 1:power])
            ]
        else
            error("Unreachable.")
        end
    else
        error("Not implemented.")
    end
end


function arcdiag_to_deformationmap__so(
    diag::ArcDiagramUndirected,
    sp::SmashProductLie{C},
    W::LieAlgebraModule=exterior_power(base_module(sp), 2),
) where {C <: RingElem}
    @req !is_direct_sum(W) "Not permitted for direct sums."
    ind_map = basis_index_mapping(W)

    # TODO: allow for genereal ArcDiagrams
    dim_stdmod_V = base_lie_algebra(sp).n

    iso_wedge2V_to_L = Dict{Vector{Int}, Int}()
    for (i, bs) in enumerate(combinations(dim_stdmod_V, 2))
        iso_wedge2V_to_L[bs] = i
    end

    if is_exterior_power(W)
        nrows_kappa = ncols_kappa = dim(base_module(W))
    elseif is_tensor_product(W)
        nrows_kappa, ncols_kappa = dim.(base_modules(W))
    end

    kappa = zero_matrix(underlying_algebra(sp), nrows_kappa, ncols_kappa)
    for (label_index, upper_labels) in enumerate(arc_diagram_label_iterator__so(W, 1:dim_stdmod_V))

        i, j = ind_map[label_index]

        entry = arcdiag_to_deformationmap_entry__so(
            diag,
            W,
            upper_labels,
            underlying_algebra(sp),
            iso_wedge2V_to_L,
            dim_stdmod_V,
        )

        entry = _normal_form(entry, sp.rels)

        kappa[i, j] += entry
        if is_exterior_power(W)
            kappa[j, i] -= entry
        end
    end
    return kappa
end

function arcdiag_to_deformationmap_entry__so(
    diag::ArcDiagramUndirected,
    W::LieAlgebraModule{C},
    upper_labels::AbstractVector{Int},
    sp_alg::FreeAssAlgebra{C},
    iso_pair_to_L::Dict{Vector{Int}, Int},
    max_label::Int,
) where {C <: RingElem}
    entry = zero(sp_alg)

    for (upper_labels, sgn_upper_labels) in arc_diagram_label_permutations__so(W, upper_labels)
        zeroprod = false
        lower_labels = [0 for _ in 1:n_lower_vertices(diag)]
        frees = Int[]
        for v in upper_vertices(diag)
            nv = neighbor(diag, v)
            if is_upper_vertex(nv) && upper_labels[vertex_index(v)] != upper_labels[vertex_index(nv)]
                zeroprod = true
                break
            elseif is_lower_vertex(nv)
                lower_labels[vertex_index(nv)] = upper_labels[vertex_index(v)]
            end
        end
        if zeroprod
            continue
        end
        for v in lower_vertices(diag)
            nv = neighbor(diag, v)
            if is_lower_vertex(nv) && vertex_index(v) < vertex_index(nv)
                push!(frees, vertex_index(v))
            end
        end

        entry_summand = zero(sp_alg)

        # iterate over lower point labelings
        nextindex = 1
        while true
            if nextindex > length(frees)
                # begin inner
                zeroelem = false
                sign_lower_labels = 1
                basiselem = Int[]
                for k in 1:2:length(lower_labels)
                    if lower_labels[k] == lower_labels[k+1]
                        zeroelem = true
                        break
                    elseif lower_labels[k] > lower_labels[k+1]
                        sign_lower_labels *= -1
                        append!(basiselem, iso_pair_to_L[[lower_labels[k+1], lower_labels[k]]])
                    else
                        append!(basiselem, iso_pair_to_L[[lower_labels[k], lower_labels[k+1]]])
                    end
                end
                if !zeroelem
                    symm_basiselem = sp_alg(
                        fill(C(1 // factorial(length(basiselem))), factorial(length(basiselem))),
                        [ind for ind in permutations(basiselem)],
                    )
                    entry_summand += sign_lower_labels * symm_basiselem # TODO: benchmark removal of normal_form
                end
                # end inner

                nextindex -= 1
            end

            while nextindex >= 1 && lower_labels[frees[nextindex]] == max_label
                lower_labels[frees[nextindex]] = 0
                lower_labels[vertex_index(_neighbor_of_lower_vertex(diag, frees[nextindex]))] = 0
                nextindex -= 1
            end
            if nextindex == 0
                break
            end
            lower_labels[frees[nextindex]] += 1
            lower_labels[vertex_index(_neighbor_of_lower_vertex(diag, frees[nextindex]))] += 1
            if ispairgood(lower_labels, frees[nextindex]) &&
               ispairgood(lower_labels, vertex_index(_neighbor_of_lower_vertex(diag, frees[nextindex])))
                nextindex += 1
            end
        end

        entry += sgn_upper_labels * entry_summand
    end
    return entry
end

function ispairgood(labeled_diag::Vector{Int}, k::Int)
    left_k = k % 2 == 1 ? k : k - 1
    return labeled_diag[left_k] != labeled_diag[left_k+1]
end
