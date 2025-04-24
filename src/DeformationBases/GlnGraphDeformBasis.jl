const GlnGraphDeformBasisDataT = Tuple{GlnGraph, Vector{Int}, Partition{Int}}

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
    extra_data::Dict{DeformationMap{T}, Set{Tuple{Tuple{Int, Int}, GlnGraphDeformBasisDataT}}}
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
        extra_data = Dict{DeformationMap{elem_type(sp)}, Set{Tuple{Tuple{Int, Int}, GlnGraphDeformBasisDataT}}}()

        function data_iter_and_len(LieType::GL, W::LieAlgebraModule, case::Symbol, d::Int)
            parity_verts = arc_diagram_upper_points(LieType, W)
            if case == :exterior_power
                fl, Wbase, k = _is_exterior_power(W)
                @assert fl
                @assert k == 2
                n_left_verts = n_right_verts = arc_diagram_num_upper_points(LieType, Wbase)
            elseif case == :tensor_product
                fl, W_factors = _is_tensor_product(W)
                @assert fl
                @assert length(W_factors) == 2
                n_left_verts, n_right_verts = arc_diagram_num_upper_points.(Ref(LieType), W_factors)
            else
                error("Unknown case")
            end
            @assert length(parity_verts) == n_left_verts + n_right_verts

            gln_graph_iter_len = number_of_gln_graphs(n_left_verts, n_right_verts, parity_verts)
            gln_graph_iter = all_gln_graphs(n_left_verts, n_right_verts, parity_verts)
            n_edges = div(length(parity_verts), 2)

            label_part_iter_len = number_of_gln_graph_labelings_with_partitions(n_edges, d)
            label_part_iter = gln_graph_labelings_with_partitions(n_edges, d)
            data_iter = ((g, labeling, part) for g in gln_graph_iter for (labeling, part) in label_part_iter)
            len = Int(gln_graph_iter_len * label_part_iter_len)
            return data_iter, len::Int
        end

        iter1, len1 = arc_diag_based_basis_iteration(
            LieType,
            sp,
            degs,
            extra_data,
            data_iter_and_len;
            no_normalize,
        )

        if no_normalize
            return new{elem_type(sp)}(len1, iter1, extra_data, no_normalize)
        else
            iter2, len2 = filter_independent(coefficient_ring(sp), iter1)
            return new{elem_type(sp)}(len2, iter2, extra_data, no_normalize)
        end
    end
end

function Base.iterate(i::GlnGraphDeformBasis)
    return iterate(i.iter)::Union{Tuple{eltype(i), Any}, Nothing}
end

function Base.iterate(i::GlnGraphDeformBasis, s)
    return iterate(i.iter, s)::Union{Tuple{eltype(i), Any}, Nothing}
end

Base.length(basis::GlnGraphDeformBasis) = basis.len
