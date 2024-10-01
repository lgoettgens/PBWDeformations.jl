function Base.:(==)(pg1::PseudographLabelled{T}, pg2::PseudographLabelled{T}) where {T}
    return (pg1.nv, pg1.edges) == (pg2.nv, pg2.edges)
end

function Base.hash(pg::PseudographLabelled, h::UInt)
    b = 0x1f83bacf6c93c455 % UInt
    h = hash(pg.nv, h)
    h = hash(edges(pg), h)
    return xor(h, b)
end

function n_vertices(pg::PseudographLabelled)
    return pg.nv
end

function n_edges(pg::PseudographLabelled)
    return Int(length(edges(pg)))
end

function edges(pg::PseudographLabelled)
    return pg.edges
end

function n_edges(pg::PseudographLabelled, verts::MSet{Int})
    return Int(length(edges(pg, verts)))
end

function edges(pg::PseudographLabelled, verts::MSet{Int})
    return filter(e -> first(e) == verts, edges(pg))
end

function edge_labels(pg::PseudographLabelled, verts::MSet{Int})
    return [last(e) for e in filter(e -> first(e) == verts, edges(pg))]
end

function Base.sum(pg::PseudographLabelled{T}) where {T <: Number}
    return sum(last, edges(pg); init=zero(T))
end

function Base.sum(pg::PseudographLabelled{T}, verts::MSet{Int}) where {T <: Number}
    return sum(edge_labels(pg, verts); init=zero(T))
end

function all_pseudographs(nv::Int, degree::Int, sumtotal::Int; upto_iso::Bool=false)
    @req nv == 2 "Only implemented for 2 vertices"

    result = PseudographLabelled{Int}[]
    for nloops in 0:div(degree, 2)
        nedges = degree - 2 * nloops
        for sumedges in 0:sumtotal
            if nedges == 0 && sumedges > 0
                continue
            end
            if nloops == 0 && sumedges < sumtotal
                continue
            end
            for sumloop2 in 0:sumtotal-sumedges
                sumloop1 = sumtotal - sumedges - sumloop2

                for edge_weights in partitions(sumedges + nedges, nedges)
                    for loop1_weights in partitions(sumloop1 + nloops, nloops)
                        for loop2_weights in partitions(sumloop2 + nloops, nloops)
                            if upto_iso && !(loop1_weights >= loop2_weights)
                                continue
                            end
                            edges = MSet(
                                [
                                    [MSet([1, 1]) => k for k in loop1_weights .- 1]
                                    [MSet([2, 2]) => k for k in loop2_weights .- 1]
                                    [MSet([1, 2]) => k for k in edge_weights .- 1]
                                ],
                            )
                            push!(result, PseudographLabelled(2, edges; regular_degree=degree, check=false))
                        end
                    end
                end
            end
        end
    end

    return result
end

function to_arcdiag(pg::PseudographLabelled{Int}, part::Partition{Int}=Partition(Int[]))
    @req n_vertices(pg) == 2 "Only implemented for 2 vertices"

    n_upper_verts = 2 * n_edges(pg)
    n_lower_verts = 2 * (sum(pg) + sum(part))

    upper_adj = zeros(Int, n_upper_verts)
    lower_adj = zeros(Int, n_lower_verts)
    addarc_uu(i, j) = (upper_adj[i] = -j; upper_adj[j] = -i)
    addarc_ll(i, j) = (lower_adj[i] = j; lower_adj[j] = i)
    addarc_ul(i, j) = (upper_adj[i] = j; lower_adj[j] = -i)
    addarc_lu(i, j) = (lower_adj[i] = -j; upper_adj[j] = i)
    i = 1
    j = div(n_upper_verts, 2) + 1
    k = 1
    for l1 in edge_labels(pg, MSet([1, 1]))
        if l1 == 0
            addarc_uu(i, i + 1)
            i += 2
        else
            addarc_ul(i, k)
            i += 1
            k += 1
            for _ in 2:l1
                addarc_ll(k, k + 1)
                k += 2
            end
            addarc_lu(k, i)
            i += 1
            k += 1
        end
    end

    for e in edge_labels(pg, MSet([1, 2]))
        if e == 0
            addarc_uu(i, j)
            i += 1
            j += 1
        else
            addarc_ul(i, k)
            i += 1
            k += 1
            for _ in 2:e
                addarc_ll(k, k + 1)
                k += 2
            end
            addarc_lu(k, j)
            j += 1
            k += 1
        end
    end

    for l2 in edge_labels(pg, MSet([2, 2]))
        if l2 == 0
            addarc_uu(j, j + 1)
            j += 2
        else
            addarc_ul(j, k)
            j += 1
            k += 1
            for _ in 2:l2
                addarc_ll(k, k + 1)
                k += 2
            end
            addarc_lu(k, j)
            j += 1
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

    return arc_diagram(Undirected, upper_adj, lower_adj)
end
