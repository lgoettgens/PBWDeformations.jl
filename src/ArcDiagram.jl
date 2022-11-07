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

    if VERSION >= v"1.7"
        function ArcDiagram(sUpper::AbstractString, sLower::AbstractString)
            sUpper = strip(sUpper)
            sLower = strip(sLower)
            nUpper = length(sUpper)
            nLower = length(sLower)

            str = sUpper * sLower
            adj = [0 for i in 1:nUpper+nLower]
            symbols = unique(str)
            for s in symbols
                count(s, str) == 2 || throw(ArgumentError("Symbol $s does not appear exactly twice."))
            end
            for s in symbols
                i, j = findall(s, str)
                adj[i] = j
                adj[j] = i
            end
            return new(nUpper, nLower, adj)
        end
    else
        function ArcDiagram(sUpper::AbstractString, sLower::AbstractString)
            sUpper = strip(sUpper)
            sLower = strip(sLower)
            nUpper = length(sUpper)
            nLower = length(sLower)

            str = sUpper * sLower
            adj = [0 for i in 1:nUpper+nLower]
            symbols = unique(str)
            for s in symbols
                count(string(s), str) == 2 || throw(ArgumentError("Symbol $s does not appear exactly twice."))
            end
            for s in symbols
                i, j = findall(string(s), str)
                adj[i] = j
                adj[j] = i
            end
            return new(nUpper, nLower, adj)
        end
    end

    function ArcDiagram(s::AbstractString)
        sUpper, sLower = split(s, ",")
        return ArcDiagram(sUpper, sLower)
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
    symbols = join(vcat('A':'Z', 'a':'z', '0':'9'))
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

function is_crossing_free(a::ArcDiagram; part=:everything::Symbol)
    if part == :everything
        for i in 1:a.nUpper+a.nLower, j in 1:a.nUpper+a.nLower
            if i == j
                continue
            elseif i >= a.adjacency[i]
                continue
            elseif j >= a.adjacency[j]
                continue
            end
            i, j = min(i, j), max(i, j)
            # now i < j
            if a.adjacency[i] <= a.nUpper
                if (j < a.adjacency[i]) != (a.adjacency[j] < a.adjacency[i])
                    return false
                end
            elseif i > a.nUpper
                if (j < a.adjacency[i]) != (a.adjacency[j] < a.adjacency[i])
                    return false
                end
            elseif a.adjacency[j] <= a.nUpper
                continue
            elseif a.nUpper < j
                if j < a.adjacency[i] < a.adjacency[j]
                    return false
                end
            elseif a.adjacency[i] > a.adjacency[j]
                return false
            end
        end
        return true
    elseif part == :upper
        for i in 1:a.nUpper, j in 1:a.nUpper
            if i == j
                continue
            elseif i >= a.adjacency[i]
                continue
            elseif j >= a.adjacency[j]
                continue
            elseif a.adjacency[i] > a.nUpper
                continue
            elseif a.adjacency[j] > a.nUpper
                continue
            end
            i, j = min(i, j), max(i, j)
            # now i < j
            if (j < a.adjacency[i]) != (a.adjacency[j] < a.adjacency[i])
                return false
            end
        end
        return true
    elseif part == :lower
        for i in 1:a.nLower, j in 1:a.nLower
            i += a.nUpper
            j += a.nUpper
            if i == j
                continue
            elseif i >= a.adjacency[i]
                continue
            elseif j >= a.adjacency[j]
                continue
            end
            i, j = min(i, j), max(i, j)
            # now i < j
            if (j < a.adjacency[i]) != (a.adjacency[j] < a.adjacency[i])
                return false
            end
        end
        return true
    else
        error("Unknown part")
    end
end


struct ArcDiagramIterator
    iter
    len::Int
end

function Base.iterate(i::ArcDiagramIterator)
    return iterate(i.iter)
end

function Base.iterate(i::ArcDiagramIterator, s)
    return iterate(i.iter, s)
end

Base.length(i::ArcDiagramIterator) = i.len

Base.eltype(::Type{ArcDiagramIterator}) = ArcDiagram


function all_arc_diagrams(nUpper::Int, nLower::Int; indep_sets::AbstractVector{<:AbstractVector{Int}}=Vector{Int}[])
    n = nUpper + nLower
    iter, len = iter_possible_adjacencies(nUpper, nLower, indep_sets, [0 for i in 1:n])
    return ArcDiagramIterator(iter, len)
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
        return [ArcDiagram(nUpper, nLower, partial; skipchecks=true)], 1
    end
    rel_indep_sets = filter(is -> i in is, indep_sets)
    poss_adjs = setdiff(setdiff(findall(==(0), partial), i), rel_indep_sets...)
    choices = Iterators.map(poss_adjs) do j
        partial2 = deepcopy(partial)
        partial2[i] = j
        partial2[j] = i
        return iter_possible_adjacencies(nUpper, nLower, indep_sets, partial2)
    end
    return Iterators.flatten(Iterators.map(c -> c[1], choices)), sum(c -> c[2], choices; init=0)
end
