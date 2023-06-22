const ArcDiagramNode = Tuple{Symbol, Int}

function is_upper_node(p::ArcDiagramNode)
    return p[1] == :upper
end

function is_lower_node(p::ArcDiagramNode)
    return p[1] == :lower
end

function node_index(p::ArcDiagramNode)
    return p[2]
end

struct ArcDiagram
    num_upper_nodes::Int
    num_lower_nodes::Int
    upper_neighbors::Vector{ArcDiagramNode}
    lower_neighbors::Vector{ArcDiagramNode}

    function ArcDiagram(
        num_upper_nodes::Int,
        num_lower_nodes::Int,
        upper_neighbors::Vector{ArcDiagramNode},
        lower_neighbors::Vector{ArcDiagramNode};
        check::Bool=true,
    )
        if check
            @req length(upper_neighbors) == num_upper_nodes "Upper nodes' neighbors list has wrong length."
            @req length(lower_neighbors) == num_lower_nodes "Lower nodes' neighbors list has wrong length."
            for i in 1:num_upper_nodes
                neigh = upper_neighbors[i]
                if is_upper_node(neigh)
                    @req 1 <= node_index(neigh) <= num_upper_nodes "Out of bounds adjacency."
                    @req node_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = upper_neighbors[node_index(neigh)]
                    @req is_upper_node(neigh_neigh) && node_index(neigh_neigh) == i "Adjacency not symmetric."
                elseif is_lower_node(neigh)
                    @req 1 <= node_index(neigh) <= num_lower_nodes "Out of bounds adjacency."
                    neigh_neigh = lower_neighbors[node_index(neigh)]
                    @req is_upper_node(neigh_neigh) && node_index(neigh_neigh) == i "Adjacency not symmetric."
                else
                    @req false "Invalid neighbor type."
                end
            end
            for i in 1:num_lower_nodes
                neigh = lower_neighbors[i]
                if is_upper_node(neigh)
                    @req 1 <= node_index(neigh) <= num_upper_nodes "Out of bounds adjacency."
                    neigh_neigh = upper_neighbors[node_index(neigh)]
                    @req is_lower_node(neigh_neigh) && node_index(neigh_neigh) == i "Adjacency not symmetric."
                elseif is_lower_node(neigh)
                    @req 1 <= node_index(neigh) <= num_lower_nodes "Out of bounds adjacency."
                    @req node_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = lower_neighbors[node_index(neigh)]
                    @req is_lower_node(neigh_neigh) && node_index(neigh_neigh) == i "Adjacency not symmetric."
                else
                    @req false "Invalid neighbor type."
                end
            end
        end
        return new(num_upper_nodes, num_lower_nodes, upper_neighbors, lower_neighbors)
    end

    function ArcDiagram(
        num_upper_nodes::Int,
        num_lower_nodes::Int,
        upper_neighbor_inds::Vector{Int},
        lower_neighbor_inds::Vector{Int};
        check::Bool=true,
    )
        upper_neighbors = [i < 0 ? (:upper, -i) : (:lower, i) for i in upper_neighbor_inds]
        lower_neighbors = [i < 0 ? (:upper, -i) : (:lower, i) for i in lower_neighbor_inds]
        return ArcDiagram(num_upper_nodes, num_lower_nodes, upper_neighbors, lower_neighbors; check)
    end

    function ArcDiagram(upper_neighbor_inds::Vector{Int}, lower_neighbor_inds::Vector{Int}; check::Bool=true)
        num_upper_nodes = length(upper_neighbor_inds)
        num_lower_nodes = length(lower_neighbor_inds)
        return ArcDiagram(num_upper_nodes, num_lower_nodes, upper_neighbor_inds, lower_neighbor_inds; check)
    end

    function ArcDiagram(upper_neighbors::Vector{Int}, lower_neighbors::Vector{Int}; check::Bool=true)
        num_upper_nodes = length(upper_neighbors)
        num_lower_nodes = length(lower_neighbors)
        return ArcDiagram(num_upper_nodes, num_lower_nodes, upper_neighbors, lower_neighbors; check)
    end

    function ArcDiagram(upper::AbstractString, lower::AbstractString)
        upper = strip(upper)
        lower = strip(lower)
        num_upper_nodes = length(upper)
        num_lower_nodes = length(lower)

        str = upper * lower
        symbols = unique(str)
        for s in symbols
            @static if VERSION >= v"1.7"
                @req count(s, str) == 2 "Symbol $s does not appear exactly twice."
            else
                @req count(string(s), str) == 2 "Symbol $s does not appear exactly twice."
            end
        end
        upper_neighbors = [(:none, 0) for _ in 1:num_upper_nodes]
        lower_neighbors = [(:none, 0) for _ in 1:num_lower_nodes]
        for s in symbols
            @static if VERSION >= v"1.7"
                i, j = findall(s, str)
            else
                i, j = findall(string(s), str)
            end
            if i <= num_upper_nodes
                upper_neighbors[i] = j <= num_upper_nodes ? (:upper, j) : (:lower, j - num_upper_nodes)
            else
                lower_neighbors[i-num_upper_nodes] = j <= num_upper_nodes ? (:upper, j) : (:lower, j - num_upper_nodes)
            end
            if j <= num_upper_nodes
                upper_neighbors[j] = i <= num_upper_nodes ? (:upper, i) : (:lower, i - num_upper_nodes)
            else
                lower_neighbors[j-num_upper_nodes] = i <= num_upper_nodes ? (:upper, i) : (:lower, i - num_upper_nodes)
            end
        end
        return new(num_upper_nodes, num_lower_nodes, upper_neighbors, lower_neighbors)
    end

    function ArcDiagram(s::AbstractString)
        upper, lower = split(s, ",")
        return ArcDiagram(upper, lower)
    end
end

function Base.:(==)(a1::ArcDiagram, a2::ArcDiagram)
    return (a1.num_upper_nodes, a1.num_lower_nodes, a1.upper_neighbors, a1.lower_neighbors) ==
           (a2.num_upper_nodes, a2.num_lower_nodes, a2.upper_neighbors, a2.lower_neighbors)
end

function Base.hash(a::ArcDiagram, h::UInt)
    b = 0x7e9086d4b4dae57d % UInt
    h = hash(a.num_upper_nodes, h)
    h = hash(a.num_lower_nodes, h)
    h = hash(a.upper_neighbors, h)
    h = hash(a.lower_neighbors, h)
    return xor(h, b)
end

function Base.show(io::IO, a::ArcDiagram)
    show(IOContext(io, (:compact => true)), MIME"text/plain"(), a)
end

function Base.show(io::IO, ::MIME"text/plain", a::ArcDiagram)
    pair_ids = map(i -> min(i, a.adjacency[i]), 1:a.num_upper_nodes+a.num_lower_nodes)
    symbols = join(vcat('A':'Z', 'a':'z', '0':'9'))
    max_id = maximum(pair_ids)
    if max_id > length(symbols)
        print(io, "too large to print")
        return
    end
    print(io, symbols[pair_ids[1:a.num_upper_nodes]])
    if get(io, :compact, false)
        print(io, ",")
    else
        print(io, "\n")
    end
    print(io, symbols[pair_ids[a.num_upper_nodes+1:a.num_upper_nodes+a.num_lower_nodes]])
end

function is_crossing_free(a::ArcDiagram; part=:everything::Symbol)
    return false # TODO: fixme
    if part == :everything
        for i in 1:a.num_upper_nodes+a.num_lower_nodes, j in 1:a.num_upper_nodes+a.num_lower_nodes
            if i == j
                continue
            elseif i >= a.adjacency[i]
                continue
            elseif j >= a.adjacency[j]
                continue
            end
            i, j = min(i, j), max(i, j)
            # now i < j
            if a.adjacency[i] <= a.num_upper_nodes
                if (j < a.adjacency[i]) != (a.adjacency[j] < a.adjacency[i])
                    return false
                end
            elseif i > a.num_upper_nodes
                if (j < a.adjacency[i]) != (a.adjacency[j] < a.adjacency[i])
                    return false
                end
            elseif a.adjacency[j] <= a.num_upper_nodes
                continue
            elseif a.num_upper_nodes < j
                if j < a.adjacency[i] < a.adjacency[j]
                    return false
                end
            elseif a.adjacency[i] > a.adjacency[j]
                return false
            end
        end
        return true
    elseif part == :upper
        for i in 1:a.num_upper_nodes, j in 1:a.num_upper_nodes
            if i == j
                continue
            elseif i >= a.adjacency[i]
                continue
            elseif j >= a.adjacency[j]
                continue
            elseif a.adjacency[i] > a.num_upper_nodes
                continue
            elseif a.adjacency[j] > a.num_upper_nodes
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
        for i in 1:a.num_lower_nodes, j in 1:a.num_lower_nodes
            i += a.num_upper_nodes
            j += a.num_upper_nodes
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


function all_arc_diagrams(
    num_upper_nodes::Int,
    num_lower_nodes::Int;
    indep_sets::AbstractVector{<:AbstractVector{Int}}=Vector{Int}[],
)
    return ArcDiagram[] # TODO: fixme

    n = num_upper_nodes + num_lower_nodes
    iter, len = iter_possible_adjacencies(num_upper_nodes, num_lower_nodes, indep_sets, [0 for i in 1:n])
    return ArcDiagramIterator(iter, len)
end

function iter_possible_adjacencies(
    num_upper_nodes::Int,
    num_lower_nodes::Int,
    indep_sets::AbstractVector{<:AbstractVector{Int}},
    partial::Vector{Int},
)
    n = num_upper_nodes + num_lower_nodes
    i = findfirst(==(0), partial)
    if i === nothing
        return [ArcDiagram(num_upper_nodes, num_lower_nodes, partial; skipchecks=true)], 1
    end
    rel_indep_sets = filter(is -> i in is, indep_sets)
    poss_adjs = setdiff(setdiff(findall(==(0), partial), i), rel_indep_sets...)
    choices = Iterators.map(poss_adjs) do j
        partial2 = deepcopy(partial)
        partial2[i] = j
        partial2[j] = i
        return iter_possible_adjacencies(num_upper_nodes, num_lower_nodes, indep_sets, partial2)
    end
    return Iterators.flatten(Iterators.map(c -> c[1], choices)), sum(c -> c[2], choices; init=0)
end
