struct ArcDiagram
    nUpper::Int
    nLower::Int
    adjacency::Vector{Int}

    function ArcDiagram(nUpper::Int, nLower::Int, adjacency::Vector{Int}; skipchecks::Bool=false)
        if !skipchecks
            n = nUpper + nLower
            length(adjacency) == n || throw(ArgumentError("Adjacency list has wrong length."))
            all(i -> 1 <= adjacency[i] <= n, 1:n) || throw(ArgumentError("Out of bounds adjacency."))
            all(i -> i != adjacency[i], 1:n) || throw(ArgumentError("No self-loops allowed."))
            all(i -> i == adjacency[adjacency[i]], 1:n) || throw(ArgumentError("Adjacency not symmetric."))
        end
        return new(nUpper, nLower, adjacency)
    end
end

function Base.:(==)(a1::ArcDiagram, a2::ArcDiagram)
    return (a1.nUpper, a1.nLower, a1.adjacency) == (a2.nUpper, a2.nLower, a2.adjacency)
end

function Base.show(io::IO, a::ArcDiagram)
    show(IOContext(io, (:compact => true)), MIME"text/plain"(), a)
end

function Base.show(io::IO, ::MIME"text/plain", a::ArcDiagram)
    pair_ids = map(i -> min(i, a.adjacency[i]), 1:a.nUpper+a.nLower)
    symbols = join(vcat('0':'9', 'A':'Z', 'a':'z'))
    max_id = maximum(pair_ids)
    if max_id > length(symbols)
        print(io, "too large to print")
        return
    end
    print(io, symbols[pair_ids[1:a.nUpper]])
    if get(io, :compact, false)
        print(io, ",")
    else
        print(io, "\n")
    end
    print(io, symbols[pair_ids[a.nUpper+1:a.nUpper+a.nLower]])
end

function all_arc_diagrams(nUpper::Int, nLower::Int; indep_sets::AbstractVector{<:AbstractVector{Int}}=Vector{Int}[])
    n = nUpper + nLower
    iter, len = iter_possible_adjacencies(nUpper, nLower, indep_sets, [0 for i in 1:n])
    return iter
end

function iter_possible_adjacencies(
    nUpper::Int,
    nLower::Int,
    indep_sets::AbstractVector{<:AbstractVector{Int}},
    partial::Vector{Int},
)
    n = nUpper + nLower
    i = findfirst(==(0), partial)
    if i === nothing
        return [ArcDiagram(nUpper, nLower, partial)], 1
    end
    rel_indep_sets = filter(is -> in(i, is), indep_sets)
    poss_adjs = setdiff(setdiff(findall(==(0), partial), i), rel_indep_sets...)
    total_iter = []
    total_len = 0
    choices = Iterators.map(poss_adjs) do j
        partial2 = deepcopy(partial)
        partial2[i] = j
        partial2[j] = i
        return iter_possible_adjacencies(nUpper, nLower, indep_sets, partial2)
    end
    return Iterators.flatten(Iterators.map(c -> c[1], choices)), sum(c -> c[2], choices; init=0)
end
