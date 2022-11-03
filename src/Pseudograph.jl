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
