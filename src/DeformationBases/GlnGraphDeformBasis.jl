"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by a GlnGraph with certain properties,
which gets transformed to an arc diagram and then handled as
in [`ArcDiagDeformBasis`](@ref).
This process is a generalization of [FM22](@cite).
"""
struct GlnGraphDeformBasis{T <: SmashProductLieElem} <: DeformBasis{T}
    len::Int
    iter
    extra_data::Dict{DeformationMap{T}, Set{Tuple{GlnGraph, Vector{Int}, Partition{Int}}}}
    no_normalize::Bool

    function GlnGraphDeformBasis(
        sp::SmashProductLie{C, LieC, LieT},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        LieType = Val(get_attribute(base_lie_algebra(sp), :type, nothing))
        @req LieType isa GL "Only works for gl_n."
        return GlnGraphDeformBasis(LieType, sp, degs; no_normalize)
    end

    function GlnGraphDeformBasis(
        LieType::GL,
        sp::SmashProductLie{C, LieC, LieT},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        V = base_module(sp)
        W = exterior_power_obj(V, 2)
        case = :exterior_power

        parity_verts = arc_diagram_upper_points(LieType, W)
        n_left_verts = n_right_verts = div(length(parity_verts), 2) # FIXME

        gln_graph_iter_len = number_of_gln_graphs(n_left_verts, n_right_verts, parity_verts)
        gln_graph_iter = all_gln_graphs(n_left_verts, n_right_verts, parity_verts)
        n_edges = div(length(parity_verts), 2)

        extra_data = Dict{DeformationMap{elem_type(sp)}, Set{Tuple{GlnGraph, Vector{Int}, Partition{Int}}}}()
        lens = []
        iters = []
        debug_counter = 0
        for d in degs
            label_part_iter_len = number_of_gln_graph_labelings_with_partitions(n_edges, d)
            label_part_iter = gln_graph_labelings_with_partitions(n_edges, d)
            iter = (
                begin
                    @vprintln :PBWDeformations 2 "Basis generation deg $(lpad(d, maximum(ndigits, degs))), $(lpad(floor(Int, 100*(debug_counter = (debug_counter % len) + 1) / len), 3))%, $(lpad(debug_counter, ndigits(len)))/$(len)"
                    diag = arc_diagram(g, labeling, part)
                    basis_elem = arcdiag_to_deformationmap(LieType, diag, sp, W, case)
                    if !no_normalize
                        basis_elem = normalize(basis_elem)
                    end
                    if haskey(extra_data, basis_elem)
                        push!(extra_data[basis_elem], (g, labeling, part))
                    else
                        extra_data[basis_elem] = Set([(g, labeling, part)])
                    end
                    basis_elem
                end for g in gln_graph_iter for (labeling, part) in label_part_iter
            )
            push!(lens, gln_graph_iter_len * label_part_iter_len)
            push!(iters, iter)
        end
        len = sum(lens)
        iter = Iterators.flatten(iters)
        if !no_normalize
            iter = unique(Iterators.filter(b -> !iszero(b), iter))
            collected = Vector{DeformationMap{elem_type(sp)}}(collect(iter))::Vector{DeformationMap{elem_type(sp)}}
            _, rels = is_linearly_independent_with_relations(coefficient_ring(sp), collected)
            inds = [findlast(!iszero, vec(rels[i, :]))::Int for i in 1:nrows(rels)]
            deleteat!(collected, inds)
            return new{elem_type(sp)}(length(collected), collected, extra_data, no_normalize)
        end
        return new{elem_type(sp)}(len, iter, extra_data, no_normalize)
    end
end

function Base.iterate(i::GlnGraphDeformBasis)
    return iterate(i.iter)::Union{Tuple{eltype(i), Any}, Nothing}
end

function Base.iterate(i::GlnGraphDeformBasis, s)
    return iterate(i.iter, s)::Union{Tuple{eltype(i), Any}, Nothing}
end

Base.length(basis::GlnGraphDeformBasis) = basis.len
