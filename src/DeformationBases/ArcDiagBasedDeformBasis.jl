struct ArcDiagBasedDeformBasis{DataT, T <: SmashProductLieElem} <: DeformBasis{T}
    len::Int
    iter
    extra_data::Dict{DeformationMap{T}, Set{Tuple{Tuple{Int, Int}, DataT}}}
    no_normalize::Bool

    function ArcDiagBasedDeformBasis{DataT}(
        sp::SmashProductLie,
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {DataT}
        LieType = Val(get_attribute(base_lie_algebra(sp), :type, nothing)::Union{Nothing, Symbol})
        @req !isnothing(LieType) "The type of Lie algebra cannot be deduced."
        return ArcDiagBasedDeformBasis{DataT}(LieType, sp, degs; no_normalize)
    end

    function ArcDiagBasedDeformBasis{DataT}(
        LieType::Union{SO, GL},
        sp::SmashProductLie,
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {DataT}
        check_input(ArcDiagBasedDeformBasis{DataT}, LieType, sp)
        @req coefficient_ring(sp) === coefficient_ring(base_lie_algebra(sp)) "Deformation bases don't support extension of the coefficient ring of the smash product."

        extra_data = Dict{DeformationMap{elem_type(sp)}, Set{Tuple{Tuple{Int, Int}, DataT}}}()

        iter1, len1 =
            arc_diag_based_basis_iteration(ArcDiagBasedDeformBasis{DataT}, LieType, sp, degs, extra_data; no_normalize)

        if no_normalize
            return new{DataT, elem_type(sp)}(len1, iter1, extra_data, no_normalize)
        else
            iter2, len2 = filter_independent(coefficient_ring(sp), iter1)
            return new{DataT, elem_type(sp)}(len2, iter2, extra_data, no_normalize)
        end
    end
end

function Base.iterate(i::ArcDiagBasedDeformBasis)
    return iterate(i.iter)::Union{Tuple{eltype(i), Any}, Nothing}
end

function Base.iterate(i::ArcDiagBasedDeformBasis, s)
    return iterate(i.iter, s)::Union{Tuple{eltype(i), Any}, Nothing}
end

Base.length(basis::ArcDiagBasedDeformBasis) = basis.len


# fallbacks
function should_use_data(
    ::Type{<:ArcDiagBasedDeformBasis},
    ::Union{SO, GL},
    ::Any,
    ::SmashProductLie,
    ::LieAlgebraModule,
    ::Symbol,
    cache::Union{Dict{<:Any, Bool}, Nothing},
)
    return true
end

function should_use_data_cache_type(::Type{<:ArcDiagBasedDeformBasis})
    return Nothing
end

function should_use_arc_diagram(
    ::Type{<:ArcDiagBasedDeformBasis},
    ::Union{SO, GL},
    ::ArcDiagram,
    ::SmashProductLie,
    ::LieAlgebraModule,
    ::Symbol,
)
    return true
end


function arc_diag_based_basis_iteration(
    ::Type{ArcDiagBasedDeformBasis{DataT}},
    LieType::Union{SO, GL},
    sp::SmashProductLie{C, LieC, LieT},
    degs::AbstractVector{Int},
    extra_data::Dict{DeformationMap{SmashProductLieElem{C, LieC, LieT}}, Set{Tuple{Tuple{Int, Int}, DataT}}};
    no_normalize::Bool,
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}, DataT}
    V = base_module(sp)

    V_nice, h = isomorphic_module_with_simple_structure(V)
    fl, V_nice_summands = _is_direct_sum(V_nice)
    if !fl
        temp = direct_sum(V_nice)
        h = compose(h, hom(V_nice, temp, identity_matrix(coefficient_ring(temp), dim(temp))))
        V_nice_summands = [V_nice]
        V_nice = temp
    end

    n_sum_cases = div(length(V_nice_summands) * (length(V_nice_summands) + 1), 2)
    lens = Int[]
    iters = Vector{DeformationMap{elem_type(sp)}}[]
    for d in degs
        sum_case = 0
        for (i_l, V_nice_summand_i_l) in enumerate(V_nice_summands),
            (i_r, V_nice_summand_i_r) in enumerate(V_nice_summands)

            if i_l > i_r
                continue
            end

            sum_case += 1

            proj_to_summand_l = compose(h, canonical_projection(V_nice, i_l))
            proj_to_summand_r = compose(h, canonical_projection(V_nice, i_r))

            if i_l == i_r
                W = exterior_power_obj(V_nice_summand_i_l, 2)
                case = :exterior_power
            else
                W = tensor_product(V_nice_summand_i_l, V_nice_summand_i_r)
                case = :tensor_product
            end

            data_iter, len = data_iter_and_len(ArcDiagBasedDeformBasis{DataT}, LieType, W, case, d)
            prog_meter = ProgressMeter.Progress(
                len;
                output=stderr,
                enabled=SHOW_PROGRESS_BARS(),
                desc="Basis generation: deg $d, case $(sum_case)/$(n_sum_cases)",
            )
            generate_showvalues(counter, data) = () -> [("iteration", (counter, data))]

            should_use_data_cache = should_use_data_cache_type(ArcDiagBasedDeformBasis{DataT})()
            @assert isnothing(should_use_data_cache) || should_use_data_cache isa Dict{<:Any, Bool}

            iter = let W=W, case=case # avoid closure capture
                data_iter |>
                enumerate |>
                Fix1(Iterators.filter, arg -> begin
                    _, data = arg
                    should_use_data(ArcDiagBasedDeformBasis{DataT}, LieType, data, sp, W, case, should_use_data_cache)
                end) |>
                Fix1(Iterators.map, arg -> begin
                    counter, data = arg
                    (counter, data, to_arc_diagram(data))
                end) |>
                Fix1(Iterators.filter, arg -> begin
                    _, _, diag = arg
                    should_use_arc_diagram(ArcDiagBasedDeformBasis{DataT}, LieType, diag, sp, W, case)
                end) |>
                Fix1(
                    Iterators.map, arg -> begin
                        counter, data, diag = arg

                        # @vprintln :PBWDeformations 2 "Basis generation deg $(lpad(d, maximum(ndigits, degs))), case $(lpad(sum_case, ndigits(n_sum_cases)))/$(n_sum_cases), $(lpad(floor(Int, 100*counter / len), 3))%, $(lpad(counter, ndigits(len)))/$(len)"
                        _basis_elem = arcdiag_to_deformationmap(LieType, diag, sp, W, case)
                        basis_elem = matrix(proj_to_summand_l) * _basis_elem * transpose(matrix(proj_to_summand_r))
                        if i_l != i_r
                            basis_elem -= transpose(basis_elem)
                        end
                        @assert is_skew_symmetric(basis_elem)

                        if !no_normalize
                            basis_elem = normalize(basis_elem)
                        end
                        if haskey(extra_data, basis_elem)
                            push!(extra_data[basis_elem], ((i_l, i_r), data))
                        else
                            extra_data[basis_elem] = Set([((i_l, i_r), data)])
                        end
                        ProgressMeter.update!(prog_meter, counter; showvalues=generate_showvalues(counter, data))
                        basis_elem
                    end,
                )
            end

            # push!(lens, len)
            # push!(iters, iter)
            collected = collect(iter)
            isnothing(should_use_data_cache) || empty!(should_use_data_cache)
            push!(lens, length(collected))
            push!(iters, collected)
            ProgressMeter.finish!(prog_meter)
        end
    end
    len = sum(lens)
    iter = Iterators.flatten(iters)

    return iter, len
end

function to_arc_diagram(diag::ArcDiagram)
    return diag
end

function to_arc_diagram(data::Tuple)
    return arc_diagram(data...)
end

function to_arc_diagram(data)
    return arc_diagram(data)
end


function arc_diagram_lower_pair_to_L(::SO, dim_stdmod_V::Int)
    # L ≅ Sᵈ ⋀² V
    iso_pair_to_L = Dict{Tuple{Int, Int}, Int}()
    for (i, bs) in enumerate(combinations(dim_stdmod_V, 2) .|> NTuple{2})
        iso_pair_to_L[bs] = i
    end
    return function (k1::Int, k2::Int)
        if k1 == k2
            return 0, 0
        elseif k1 < k2
            return iso_pair_to_L[(k1, k2)], 1
        else
            return iso_pair_to_L[(k2, k1)], -1
        end
    end
end

function arc_diagram_lower_pair_to_L(::GL, dim_stdmod_V::Int)
    # L ≅ Sᵈ (V ⊗ V*)
    iso_pair_to_L = Dict{Tuple{Int, Int}, Int}()
    for (i, bs) in enumerate(ProductIterator(1:dim_stdmod_V, 2) .|> reverse .|> NTuple{2})
        iso_pair_to_L[bs] = i
    end
    return function (k1::Int, k2::Int)
        return iso_pair_to_L[(k1, k2)], 1
    end
    return iso_pair_to_L
end

function arc_diagram_label_iterator(T::Union{SO, GL}, V::LieAlgebraModule, base_labels::AbstractVector{Int})
    if is_tensor_generator(V)
        return [[l] for l in base_labels]
    elseif ((fl, inner_mods) = _is_tensor_product(V); fl)
        return ProductIterator([
                   arc_diagram_label_iterator(T, inner_mod, base_labels) for inner_mod in reverse(inner_mods)
               ]) .|>
               reverse .|>
               Iterators.flatten .|>
               collect
    elseif ((fl, inner_mod, power) = _is_exterior_power(V); fl)
        return combinations(collect(arc_diagram_label_iterator(T, inner_mod, base_labels)), power) .|>
               Iterators.flatten .|>
               collect
    elseif ((fl, inner_mod, power) = _is_symmetric_power(V); fl)
        return multicombinations(collect(arc_diagram_label_iterator(T, inner_mod, base_labels)), power) .|>
               Iterators.flatten .|>
               collect
    elseif ((fl, inner_mod, power) = _is_tensor_power(V); fl)
        return ProductIterator(arc_diagram_label_iterator(T, inner_mod, base_labels), power) .|>
               reverse .|>
               Iterators.flatten .|>
               collect
    else
        error("Not implemented.")
    end
end


function basis_index_mapping(V::LieAlgebraModule)
    if ((fl, inner_mods) = _is_tensor_product(V); fl)
        return ProductIterator([1:dim(inner_mod) for inner_mod in reverse(inner_mods)]) .|> reverse .|> collect
    elseif ((fl, inner_mod, power) = _is_exterior_power(V); fl)
        return combinations(dim(inner_mod), power) .|> collect
    elseif ((fl, inner_mod, power) = _is_symmetric_power(V); fl)
        return multicombinations(dim(inner_mod), power) .|> collect
    elseif ((fl, inner_mod, power) = _is_tensor_power(V); fl)
        return ProductIterator(1:dim(inner_mod), power) .|> reverse .|> collect
    else
        error("Not implemented.")
    end
end

function arc_diagram_label_permutations(T::Union{SO, GL}, V::LieAlgebraModule, label::AbstractVector{Int})
    G, h = acting_group_with_sgn(V)
    @req length(label) == degree(G) "Number of labels mismatch."
    return [(permuted(label, g), sign(h(g))) for g in G]
end

function arcdiag_is_lower_pair_label_bad(::SO, labeled_diag::Vector{Int}, k::Int)
    left_k = k % 2 == 1 ? k : k - 1
    return labeled_diag[left_k] == labeled_diag[left_k+1]
end

function arcdiag_is_lower_pair_label_bad(::GL, labeled_diag::Vector{Int}, k::Int)
    return false
end

function arcdiag_to_deformationmap(
    T::GL,
    diag::ArcDiagramDirected,
    sp::SmashProductLie{C},
    W::LieAlgebraModule=exterior_power_obj(base_module(sp), 2),
    case::Symbol=:exterior_power,
) where {C <: RingElem}
    return arcdiag_to_deformationmap(T, arc_diagram(Undirected, diag), sp, W, case)
end

function arcdiag_to_deformationmap(
    T::Union{SO, GL},
    diag::ArcDiagramUndirected,
    sp::SmashProductLie{C},
    W::LieAlgebraModule=exterior_power_obj(base_module(sp), 2),
    case::Symbol=:exterior_power,
) where {C <: RingElem}
    @req !_is_direct_sum(W)[1] "Not permitted for direct sums."
    ind_map = basis_index_mapping(W)

    dim_stdmod_V = base_lie_algebra(sp).n

    iso_pair_to_L = arc_diagram_lower_pair_to_L(T, dim_stdmod_V)

    if case == :exterior_power
        fl, Wbase, k = _is_exterior_power(W)
        @assert fl
        @assert k == 2
        nrows_kappa = ncols_kappa = dim(Wbase)
    elseif case == :tensor_product
        fl, W_factors = _is_tensor_product(W)
        @assert fl
        @assert length(W_factors) == 2
        nrows_kappa, ncols_kappa = dim.(W_factors)
    else
        error("Unknown case")
    end

    kappa = zero_matrix(sp, nrows_kappa, ncols_kappa)
    for (label_index, upper_labels) in enumerate(arc_diagram_label_iterator(T, W, 1:dim_stdmod_V))

        i, j = ind_map[label_index]

        entry = arcdiag_to_deformationmap_entry(T, diag, W, upper_labels, sp, iso_pair_to_L, dim_stdmod_V)
        simplify(entry)

        kappa[i, j] += entry
        if case == :exterior_power
            kappa[j, i] -= entry
        end
    end
    return kappa
end


function arcdiag_to_deformationmap_entry(
    T::Union{SO, GL},
    diag::ArcDiagramUndirected,
    W::LieAlgebraModule{C},
    upper_labels::AbstractVector{Int},
    sp::SmashProductLie{C},
    iso_pair_to_L::Function,
    max_label::Int,
) where {C <: RingElem}
    entry = zero(sp)

    for (upper_labels, sgn_upper_labels) in arc_diagram_label_permutations(T, W, upper_labels)
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

        entry_summand = zero(sp)

        # iterate over lower point labelings
        nextindex = 1
        while true
            if nextindex > length(frees)
                # begin inner
                coeff_lower_labels = 1
                basiselem = Int[]
                for k in 1:2:length(lower_labels)
                    gen_ind, coeff = iso_pair_to_L(lower_labels[k], lower_labels[k+1])
                    coeff_lower_labels *= coeff
                    if iszero(coeff_lower_labels)
                        break
                    end
                    append!(basiselem, gen_ind)
                end
                if !iszero(coeff_lower_labels)
                    symm_basiselem = sp(
                        underlying_algebra(sp)(
                            fill(coefficient_ring(sp)(1 // factorial(length(basiselem))), factorial(length(basiselem))),
                            [ind for ind in permutations(basiselem)],
                        ),
                    ) # TODO: benchmark use of `symmetrize` here once it is implemented with mutable arithmetics
                    entry_summand += coeff_lower_labels * symm_basiselem
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

            if arcdiag_is_lower_pair_label_bad(T, lower_labels, frees[nextindex]) || arcdiag_is_lower_pair_label_bad(
                T,
                lower_labels,
                vertex_index(_neighbor_of_lower_vertex(diag, frees[nextindex])),
            )
                continue
            end

            nextindex += 1
        end

        entry += sgn_upper_labels * entry_summand
    end
    return entry
end
