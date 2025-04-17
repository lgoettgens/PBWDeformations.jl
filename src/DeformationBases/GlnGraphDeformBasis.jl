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
        extra_data = Dict{DeformationMap{elem_type(sp)}, Set{Tuple{GlnGraph, Vector{Int}, Partition{Int}}}}()

        function diag_data_iter_and_len(LieType::GL, W::LieAlgebraModule, case::Symbol, d::Int)
            parity_verts = arc_diagram_upper_points(LieType, W)
            n_left_verts = n_right_verts = div(length(parity_verts), 2) # FIXME

            gln_graph_iter_len = number_of_gln_graphs(n_left_verts, n_right_verts, parity_verts)
            gln_graph_iter = all_gln_graphs(n_left_verts, n_right_verts, parity_verts)
            n_edges = div(length(parity_verts), 2)

            label_part_iter_len = number_of_gln_graph_labelings_with_partitions(n_edges, d)
            label_part_iter = gln_graph_labelings_with_partitions(n_edges, d)
            diag_data_iter = ((arc_diagram(g, labeling, part), (g, labeling, part)) for g in gln_graph_iter for (labeling, part) in label_part_iter)
            len = Int(gln_graph_iter_len * label_part_iter_len)
            return diag_data_iter, len::Int
        end

        function should_be_used(LieType::GL, diag::ArcDiagram, data)
            true
        end

        iter, len = arc_diag_based_basis_iteration(
            LieType,
            sp,
            degs,
            extra_data,
            diag_data_iter_and_len,
            should_be_used;
            no_normalize,
        )

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
