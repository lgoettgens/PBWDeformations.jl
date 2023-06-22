const ArcDiagramVertex = Tuple{Symbol, Int}

function is_upper_vertex(p::ArcDiagramVertex)
    return p[1] == :upper
end

function is_lower_vertex(p::ArcDiagramVertex)
    return p[1] == :lower
end

function vertex_index(p::ArcDiagramVertex)
    return p[2]
end

struct ArcDiagram
    num_upper_verts::Int
    num_lower_verts::Int
    upper_neighbors::Vector{ArcDiagramVertex}
    lower_neighbors::Vector{ArcDiagramVertex}

    function ArcDiagram(
        num_upper_verts::Int,
        num_lower_verts::Int,
        upper_neighbors::Vector{ArcDiagramVertex},
        lower_neighbors::Vector{ArcDiagramVertex};
        check::Bool=true,
    )
        if check
            @req length(upper_neighbors) == num_upper_verts "Upper vertices' neighbors list has wrong length."
            @req length(lower_neighbors) == num_lower_verts "Lower vertices' neighbors list has wrong length."
            for i in 1:num_upper_verts
                neigh = upper_neighbors[i]
                if is_upper_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= num_upper_verts "Out of bounds adjacency."
                    @req vertex_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = upper_neighbors[vertex_index(neigh)]
                    @req is_upper_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                elseif is_lower_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= num_lower_verts "Out of bounds adjacency."
                    neigh_neigh = lower_neighbors[vertex_index(neigh)]
                    @req is_upper_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                else
                    @req false "Invalid neighbor type."
                end
            end
            for i in 1:num_lower_verts
                neigh = lower_neighbors[i]
                if is_upper_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= num_upper_verts "Out of bounds adjacency."
                    neigh_neigh = upper_neighbors[vertex_index(neigh)]
                    @req is_lower_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                elseif is_lower_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= num_lower_verts "Out of bounds adjacency."
                    @req vertex_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = lower_neighbors[vertex_index(neigh)]
                    @req is_lower_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                else
                    @req false "Invalid neighbor type."
                end
            end
        end
        return new(num_upper_verts, num_lower_verts, upper_neighbors, lower_neighbors)
    end

    function ArcDiagram(
        num_upper_verts::Int,
        num_lower_verts::Int,
        upper_neighbor_inds::Vector{Int},
        lower_neighbor_inds::Vector{Int};
        check::Bool=true,
    )
        upper_neighbors = [i < 0 ? (:upper, -i) : (:lower, i) for i in upper_neighbor_inds]
        lower_neighbors = [i < 0 ? (:upper, -i) : (:lower, i) for i in lower_neighbor_inds]
        return ArcDiagram(num_upper_verts, num_lower_verts, upper_neighbors, lower_neighbors; check)
    end

    function ArcDiagram(upper_neighbor_inds::Vector{Int}, lower_neighbor_inds::Vector{Int}; check::Bool=true)
        num_upper_verts = length(upper_neighbor_inds)
        num_lower_verts = length(lower_neighbor_inds)
        return ArcDiagram(num_upper_verts, num_lower_verts, upper_neighbor_inds, lower_neighbor_inds; check)
    end

    function ArcDiagram(upper_neighbors::Vector{Int}, lower_neighbors::Vector{Int}; check::Bool=true)
        num_upper_verts = length(upper_neighbors)
        num_lower_verts = length(lower_neighbors)
        return ArcDiagram(num_upper_verts, num_lower_verts, upper_neighbors, lower_neighbors; check)
    end

    function ArcDiagram(upper::AbstractString, lower::AbstractString)
        upper = strip(upper)
        lower = strip(lower)
        num_upper_verts = length(upper)
        num_lower_verts = length(lower)

        str = upper * lower
        symbols = unique(str)
        for s in symbols
            @static if VERSION >= v"1.7"
                @req count(s, str) == 2 "Symbol $s does not appear exactly twice."
            else
                @req count(string(s), str) == 2 "Symbol $s does not appear exactly twice."
            end
        end
        upper_neighbors = [(:none, 0) for _ in 1:num_upper_verts]
        lower_neighbors = [(:none, 0) for _ in 1:num_lower_verts]
        for s in symbols
            @static if VERSION >= v"1.7"
                i, j = findall(s, str)
            else
                i, j = findall(string(s), str)
            end
            if i <= num_upper_verts
                upper_neighbors[i] = j <= num_upper_verts ? (:upper, j) : (:lower, j - num_upper_verts)
            else
                lower_neighbors[i-num_upper_verts] = j <= num_upper_verts ? (:upper, j) : (:lower, j - num_upper_verts)
            end
            if j <= num_upper_verts
                upper_neighbors[j] = i <= num_upper_verts ? (:upper, i) : (:lower, i - num_upper_verts)
            else
                lower_neighbors[j-num_upper_verts] = i <= num_upper_verts ? (:upper, i) : (:lower, i - num_upper_verts)
            end
        end
        return new(num_upper_verts, num_lower_verts, upper_neighbors, lower_neighbors)
    end

    function ArcDiagram(s::AbstractString)
        upper, lower = split(s, ",")
        return ArcDiagram(upper, lower)
    end
end

function Base.:(==)(a1::ArcDiagram, a2::ArcDiagram)
    return (a1.num_upper_verts, a1.num_lower_verts, a1.upper_neighbors, a1.lower_neighbors) ==
           (a2.num_upper_verts, a2.num_lower_verts, a2.upper_neighbors, a2.lower_neighbors)
end

function Base.hash(a::ArcDiagram, h::UInt)
    b = 0x7e9086d4b4dae57d % UInt
    h = hash(a.num_upper_verts, h)
    h = hash(a.num_lower_verts, h)
    h = hash(a.upper_neighbors, h)
    h = hash(a.lower_neighbors, h)
    return xor(h, b)
end

function Base.show(io::IO, a::ArcDiagram)
    show(IOContext(io, (:compact => true)), MIME"text/plain"(), a)
end

function Base.show(io::IO, ::MIME"text/plain", a::ArcDiagram)
    symbols = join(vcat('A':'Z', 'a':'z', '0':'9'))
    if a.num_upper_verts + a.num_lower_verts > length(symbols)
        print(io, "Very large arc diagram")
        return
    end
    symb_dict = Dict(v => s for (v, s) in zip(vertices(a), symbols))
    function first_vertex(v1::ArcDiagramVertex, v2::ArcDiagramVertex)
        if is_upper_vertex(v1) && is_lower_vertex(v2)
            return v1
        elseif is_lower_vertex(v1) && is_upper_vertex(v2)
            return v2
        else
            return vertex_index(v1) <= vertex_index(v2) ? v1 : v2

        end
    end
    print(io, join(symb_dict[first_vertex(v, neighbor(a, v))] for v in upper_vertices(a)))
    if get(io, :compact, false)
        print(io, ",")
    else
        print(io, "\n")
    end
    print(io, join(symb_dict[first_vertex(v, neighbor(a, v))] for v in lower_vertices(a)))
end

function vertices(a::ArcDiagram)
    return vcat(upper_vertices(a), lower_vertices(a))
end

function upper_vertices(a::ArcDiagram)
    return ArcDiagramVertex[(:upper, i) for i in 1:a.num_upper_verts]
end

function lower_vertices(a::ArcDiagram)
    return ArcDiagramVertex[(:lower, i) for i in 1:a.num_lower_verts]
end

function neighbor(a::ArcDiagram, v::ArcDiagramVertex)
    if is_upper_vertex(v)
        @req 1 <= vertex_index(v) <= a.num_upper_verts "Out of bounds access."
        return a.upper_neighbors[vertex_index(v)]
    elseif is_lower_vertex(v)
        @req 1 <= vertex_index(v) <= a.num_lower_verts "Out of bounds access."
        return a.lower_neighbors[vertex_index(v)]
    else
        @error "Invalid vertex."
    end
end

function neighbors(a::ArcDiagram, v::ArcDiagramVertex)
    return [neighbor(a, v)]
end

function is_crossing_free(a::ArcDiagram; part=:everything::Symbol)
    return false # TODO: fixme
    if part == :everything
        for i in 1:a.num_upper_verts+a.num_lower_verts, j in 1:a.num_upper_verts+a.num_lower_verts
            if i == j
                continue
            elseif i >= a.adjacency[i]
                continue
            elseif j >= a.adjacency[j]
                continue
            end
            i, j = min(i, j), max(i, j)
            # now i < j
            if a.adjacency[i] <= a.num_upper_verts
                if (j < a.adjacency[i]) != (a.adjacency[j] < a.adjacency[i])
                    return false
                end
            elseif i > a.num_upper_verts
                if (j < a.adjacency[i]) != (a.adjacency[j] < a.adjacency[i])
                    return false
                end
            elseif a.adjacency[j] <= a.num_upper_verts
                continue
            elseif a.num_upper_verts < j
                if j < a.adjacency[i] < a.adjacency[j]
                    return false
                end
            elseif a.adjacency[i] > a.adjacency[j]
                return false
            end
        end
        return true
    elseif part == :upper
        for i in 1:a.num_upper_verts, j in 1:a.num_upper_verts
            if i == j
                continue
            elseif i >= a.adjacency[i]
                continue
            elseif j >= a.adjacency[j]
                continue
            elseif a.adjacency[i] > a.num_upper_verts
                continue
            elseif a.adjacency[j] > a.num_upper_verts
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
        for i in 1:a.num_lower_verts, j in 1:a.num_lower_verts
            i += a.num_upper_verts
            j += a.num_upper_verts
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
    num_upper_verts::Int,
    num_lower_verts::Int;
    indep_sets::AbstractVector{<:AbstractVector{Int}}=Vector{Int}[],
)
    iter, len = iter_possible_adjacencies(
        num_upper_verts,
        num_lower_verts,
        [0 for _ in 1:num_upper_verts],
        [0 for _ in 1:num_lower_verts],
        indep_sets,
    )
    return ArcDiagramIterator(iter, len)
end

function iter_possible_adjacencies(
    num_upper_verts::Int,
    num_lower_verts::Int,
    partial_upper::Vector{Int},
    partial_lower::Vector{Int},
    indep_sets::AbstractVector{<:AbstractVector{Int}},
)
    i = findfirst(==(0), partial_upper)
    if !isnothing(i)
        i = -i
        relevant_indep_sets = filter(is -> i in is, indep_sets)
        poss_upper_adjs = setdiff(setdiff(map(j -> -j, findall(==(0), partial_upper)), i), relevant_indep_sets...)
        poss_lower_adjs = setdiff(findall(==(0), partial_lower), relevant_indep_sets...)
        choices = Iterators.map([poss_upper_adjs; poss_lower_adjs]) do j
            partial_upper2 = deepcopy(partial_upper)
            partial_lower2 = deepcopy(partial_lower)
            partial_upper2[-i] = j
            if j < 0
                partial_upper2[-j] = i
            else
                partial_lower2[j] = i
            end
            iter_possible_adjacencies(num_upper_verts, num_lower_verts, partial_upper2, partial_lower2, indep_sets)
        end
        return Iterators.flatten(Iterators.map(c -> c[1], choices)), sum(c -> c[2], choices; init=0)
    else
        i = findfirst(==(0), partial_lower)
        if !isnothing(i)
            relevant_indep_sets = filter(is -> i in is, indep_sets)
            poss_lower_adjs = setdiff(setdiff(findall(==(0), partial_lower), i), relevant_indep_sets...)
            choices = Iterators.map(poss_lower_adjs) do j
                partial_lower2 = deepcopy(partial_lower)
                partial_lower2[i] = j
                partial_lower2[j] = i
                iter_possible_adjacencies(num_upper_verts, num_lower_verts, partial_upper, partial_lower2, indep_sets)
            end
            return Iterators.flatten(Iterators.map(c -> c[1], choices)), sum(c -> c[2], choices; init=0)
        else
            return [ArcDiagram(num_upper_verts, num_lower_verts, partial_upper, partial_lower)], 1
        end
    end
end
