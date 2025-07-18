function n_vertices(g::GlnGraph)
    return g.n_left_verts + g.n_right_verts
end

function n_edges(g::GlnGraph)
    return length(g.edges)
end

function Base.:(==)(g1::GlnGraph, g2::GlnGraph)
    g1.n_left_verts == g2.n_left_verts || return false
    g1.n_right_verts == g2.n_right_verts || return false
    g1.parity_verts == g2.parity_verts || return false
    return g1.edges == g2.edges
end

function Base.hash(g::GlnGraph, h::UInt)
    h = hash(g.n_left_verts, h)
    h = hash(g.n_right_verts, h)
    h = hash(g.parity_verts, h)
    h = hash(g.edges, h)
    return h
end

function _lt(g1::GlnGraph, g2::GlnGraph)
    @req g1.n_left_verts == g2.n_left_verts && g1.n_right_verts == g2.n_right_verts "number of vertices mismatch"
    @req g1.parity_verts == g2.parity_verts "parity mismatch"
    return g1.edges < g2.edges
end

################################################################################
#
# Group actions
#
################################################################################

function gset_by_type(G::PermGroup, Omega, ::Type{GlnGraph}; closed::Bool = false)
    return GSetByElements(G, ^, Omega; closed = closed, check = false)
end

function Base.:^(g::GlnGraph, p::PermGroupElem)
    @req length(g.parity_verts) == degree(p) "permutation degree mismatch"
    @req g.parity_verts == permuted(g.parity_verts, p) "parity is not preserved"
    return GlnGraph(g.n_left_verts, g.n_right_verts, g.parity_verts, on_tuples.(g.edges, p); check=false, sort=true)
end

################################################################################
#
# Iterators
#
################################################################################


function all_gln_graphs(n_left_verts::Int, n_right_verts::Int, parity_verts::Vector{Bool})
    @req n_left_verts + n_right_verts == length(parity_verts) "number of vertices mismatch"

    result = GlnGraph[]
    !iszero(parity_diff(parity_verts)) && return result # no graph if parity is not balanced

    out_indices = findall(parity_verts)
    in_indices = findall(!, parity_verts)

    return (
        begin
            edges = collect(zip(out_indices, in_indices_permuted))
            GlnGraph(n_left_verts, n_right_verts, parity_verts, edges; check=false, sort=false)
        end for in_indices_permuted in permutations(in_indices)
    )
end

function number_of_gln_graphs(n_left_verts::Int, n_right_verts::Int, parity_verts::Vector{Bool})
    @req n_left_verts + n_right_verts == length(parity_verts) "number of vertices mismatch"
    !iszero(parity_diff(parity_verts)) && return 0 # no graph if parity is not balanced

    return factorial(div(length(parity_verts), 2))
end


function gln_graph_labelings_with_partitions(n_edges::Int, total_weight::Int)
    return (
        (Vector{Int}(wcomp), part) for edge_weights in 0:total_weight for part in Iterators.reverse(collect(partitions(total_weight - edge_weights))) for
        wcomp in weak_compositions(edge_weights, n_edges)
    )
end

function number_of_gln_graph_labelings_with_partitions(n_edges::Int, total_weight::Int)
    return sum(
        part_weight -> n_weak_compositions(total_weight - part_weight, n_edges) * n_partitions(part_weight),
        0:total_weight;
        init=zero(ZZ),
    )
end


function arc_diagram(g::GlnGraph, labeling::Vector{Int}, part::Partition{Int}=Partition(Int[]))
    @req n_edges(g) == length(labeling) "number of edge labels mismatch"

    n_upper_verts = 2 * n_edges(g)
    n_lower_verts = 2 * (sum(labeling) + sum(part))

    parity_upper_verts = g.parity_verts
    @assert length(parity_upper_verts) == n_upper_verts
    parity_lower_verts = fill(false, n_lower_verts)
    parity_lower_verts[1:2:end] .= true

    upper_adj = zeros(Int, n_upper_verts)
    lower_adj = zeros(Int, n_lower_verts)
    addarc_uu(i, j) = (upper_adj[i] = -j; upper_adj[j] = -i)
    addarc_ll(i, j) = (lower_adj[i] = j; lower_adj[j] = i)
    addarc_ul(i, j) = (upper_adj[i] = j; lower_adj[j] = -i)
    addarc_lu(i, j) = (lower_adj[i] = -j; upper_adj[j] = i)
    k = 1
    for ((v, w), label) in zip(g.edges, labeling)
        if label == 0
            addarc_uu(v, w)
        else
            addarc_ul(v, k)
            k += 1
            for _ in 2:label
                addarc_ll(k, k + 1)
                k += 2
            end
            addarc_lu(k, w)
            k += 1
        end
    end

    for p in part
        start = k
        k += 1
        for _ in 2:p
            addarc_ll(k, k + 1)
            k += 2
        end
        addarc_ll(k, start)
        k += 1
    end

    return arc_diagram(Directed, parity_upper_verts, parity_lower_verts, upper_adj, lower_adj)
end
