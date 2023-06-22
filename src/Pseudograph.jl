struct Pseudograph2
    loops1::Vector{Int}
    loops2::Vector{Int}
    edges::Vector{Int}

    function Pseudograph2(
        loops1::AbstractVector{Int},
        loops2::AbstractVector{Int},
        edges::AbstractVector{Int},
        check::Bool=true,
    )
        if check
            @req all(x -> x >= 0, loops1) "all labels in loops1 must be non-negative"
            @req all(x -> x >= 0, loops2) "all labels in loops2 must be non-negative"
            @req all(x -> x >= 0, edges) "all labels in edges must be non-negative"
            @req length(loops1) == length(loops2) "both vertices must have the same degree"
        end
        return new(loops1, loops2, edges)
    end
end

function Base.:(==)(pg1::Pseudograph2, pg2::Pseudograph2)
    return (pg1.loops1, pg1.loops2, pg1.edges) == (pg2.loops1, pg2.loops2, pg2.edges)
end

function Base.hash(pg::Pseudograph2, h::UInt)
    h = hash(pg.loops1, h)
    h = hash(pg.loops2, h)
    h = hash(pg.edges, h)
    return h
end

function nedges(pg::Pseudograph2)
    return length(pg.loops1) + length(pg.loops2) + length(pg.edges)
end

function Base.sum(pg::Pseudograph2)
    return sum(pg.loops1) + sum(pg.loops2) + sum(pg.edges)
end

function all_pseudographs(reg::Int, sumtotal::Int; upto_iso::Bool=false)
    result = Pseudograph2[]
    for nloops in 0:div(reg, 2)
        ne = reg - 2 * nloops
        for sume in 0:sumtotal
            if ne == 0 && sume > 0
                continue
            end
            if nloops == 0 && sume < sumtotal
                continue
            end
            for sumloop2 in 0:sumtotal-sume
                sumloop1 = sumtotal - sume - sumloop2

                for edges in partitions(sume + ne, ne)
                    for loops1 in partitions(sumloop1 + nloops, nloops)
                        for loops2 in partitions(sumloop2 + nloops, nloops)
                            if upto_iso && !(loops1 >= loops2)
                                continue
                            end
                            push!(result, Pseudograph2(loops1 .- 1, loops2 .- 1, edges .- 1))
                        end
                    end
                end
            end
        end
    end

    return result
end

function to_arcdiag(pg::Pseudograph2, part::Generic.Partition=Partition(Int[]))
    num_upper_verts = 2 * nedges(pg)
    num_lower_verts = 2 * (sum(pg) + part.n)

    upper_adj = zeros(Int, num_upper_verts)
    lower_adj = zeros(Int, num_lower_verts)
    addarc_uu(i, j) = (upper_adj[i] = -j; upper_adj[j] = -i)
    addarc_ll(i, j) = (lower_adj[i] = j; lower_adj[j] = i)
    addarc_ul(i, j) = (upper_adj[i] = j; lower_adj[j] = -i)
    addarc_lu(i, j) = (lower_adj[i] = -j; upper_adj[j] = i)
    i = 1
    j = div(num_upper_verts, 2) + 1
    k = 1
    for l1 in pg.loops1
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

    for e in pg.edges
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

    for l2 in pg.loops2
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

    for p in part.part
        start = k
        k += 1
        for _ in 2:p
            addarc_ll(k, k + 1)
            k += 2
        end
        addarc_ll(k, start)
        k += 1
    end

    return ArcDiagram(num_upper_verts, num_lower_verts, upper_adj, lower_adj)
end
