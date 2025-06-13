const GlnGraphDeformBasisDataT = Tuple{GlnGraph, Vector{Int}, Partition{Int}}

"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by a GlnGraph with certain properties,
which gets transformed to an arc diagram and then handled as
in [`ArcDiagDeformBasis`](@ref).
This process is a generalization of [FM22](@cite).
"""
const GlnGraphDeformBasis{T} = ArcDiagBasedDeformBasis{GlnGraphDeformBasisDataT, T} where {T <: SmashProductLieElem}

function check_input(
    ::Type{GlnGraphDeformBasis},
    LieType,
    sp::SmashProductLie{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req LieType isa GL "Only works for gl_n."
end

function data_iter_and_len(::Type{GlnGraphDeformBasis}, LieType::GL, W::LieAlgebraModule, case::Symbol, d::Int)
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

function should_use_data_cache_type(::Type{GlnGraphDeformBasis})
    return Dict{Tuple{GlnGraph, Vector{Int}}, Bool}
end

function should_use_data(
    ::Type{GlnGraphDeformBasis},
    LieType::GL,
    data::GlnGraphDeformBasisDataT,
    ::SmashProductLie,
    V::LieAlgebraModule,
    ::Symbol,
    cache::Union{Dict{<:Any, Bool}, Nothing},
)
    @assert cache isa should_use_data_cache_type(GlnGraphDeformBasis)

    g, labeling, part = data

    return get!(cache, (g, labeling)) do
        G, sgn = acting_group_with_sgn(V)

        g_orb = orbit(G, g)
        # graph should be minimal in its orbit
        any(g2 -> _lt(g2, g), g_orb) && return false

        # restrict action to stabilizer of the graph to examine action on the labeling
        g_stab, _ = stabilizer(g_orb)
        is_trivial(g_stab) && return true # nothing to do

        permute_edge_action = action_homomorphism(gset(g_stab, findall(g.parity_verts); closed=true))
        labeling_orb = orbit(g_stab, induced_action(permuted, permute_edge_action), labeling)

        # labeling should be maximal in its orbit
        any(>(labeling), labeling_orb) && return false

        labeling_stab, _ = stabilizer(labeling_orb)

        # signed orbit contains element twice with opposite signs
        is_trivial(image(sgn, labeling_stab)[1]) || return false

        # no criteria left to check
        return true
    end
end
