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
            all(x -> x >= 0, loops1) || error("all labels in loops1 must be non-negative")
            all(x -> x >= 0, loops2) || error("all labels in loops2 must be non-negative")
            all(x -> x >= 0, edges) || error("all labels in edges must be non-negative")
            length(loops1) == length(loops2) || error("both vertices must have the same degree")
        end
        return new(loops1, loops2, edges)
    end
end

function nedges(pg::Pseudograph2)
    return length(pg.loops1) + length(pg.loops2) + length(pg.edges)
end

function Base.sum(pg::Pseudograph2)
    return sum(pg.loops1) + sum(pg.loops2) + sum(pg.edges)
end

function all_pseudographs(reg::Int, sum::Int)
    result = Pseudograph2[]
    for nloops in 0:div(reg, 2)
        ne = reg - 2 * nloops
        for sume in 0:sum
            if ne == 0 && sume > 0
                continue
            end
            if nloops == 0 && sume < sum
                continue
            end
            for sumloop1 in 0:sum-sume
                sumloop2 = sum - sume - sumloop1

                for edges in (ne > 0 ? Combinatorics.partitions(sume + ne, ne) : [Int[]])
                    for loops1 in (nloops > 0 ? Combinatorics.partitions(sumloop1 + nloops, nloops) : [Int[]])
                        for loops2 in (nloops > 0 ? Combinatorics.partitions(sumloop2 + nloops, nloops) : [Int[]])
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
    nUpper = 2 * nedges(pg)
    nLower = 2 * (sum(pg) + part.n)

    adj = zeros(Int, nUpper + nLower)
    i = 1
    j = div(nUpper, 2) + 1
    k = nUpper + 1
    for l1 in pg.loops1
        if l1 == 0
            adj[i] = i + 1
            adj[i+1] = i
            i += 2
        else
            adj[i] = k
            adj[k] = i
            i += 1
            k += 1
            for _ in 2:l1
                adj[k] = k + 1
                adj[k+1] = k
                k += 2
            end
            adj[k] = i
            adj[i] = k
            i += 1
            k += 1
        end
    end

    for e in pg.edges
        if e == 0
            adj[i] = j
            adj[j] = i
            i += 1
            j += 1
        else
            adj[i] = k
            adj[k] = i
            i += 1
            k += 1
            for _ in 2:e
                adj[k] = k + 1
                adj[k+1] = k
                k += 2
            end
            adj[k] = j
            adj[j] = k
            j += 1
            k += 1
        end
    end

    for l2 in pg.loops2
        if l2 == 0
            adj[j] = j + 1
            adj[j+1] = j
            j += 2
        else
            adj[j] = k
            adj[k] = j
            j += 1
            k += 1
            for _ in 2:l2
                adj[k] = k + 1
                adj[k+1] = k
                k += 2
            end
            adj[k] = j
            adj[j] = k
            j += 1
            k += 1
        end
    end

    for p in part.part
        start = k
        k += 1
        for _ in 2:p
            adj[k] = k + 1
            adj[k+1] = k
            k += 2
        end
        adj[k] = start
        adj[start] = k
        k += 1
    end

    return ArcDiagram(nUpper, nLower, adj)
end
