const ArcDiagramVertex = Tuple{Symbol, Int}

abstract type ArcDiagram end


function vertices(a::ArcDiagram)
    return vcat(upper_vertices(a), lower_vertices(a))
end

function upper_vertices(a::ArcDiagram)
    return ArcDiagramVertex[upper_vertex(ArcDiagramVertex, i) for i in 1:n_upper_vertices(a)]
end

function upper_vertex(a::ArcDiagram, i::Int)
    @req 1 <= i <= n_upper_vertices(a) "Invalid index"
    return (:upper, i)::ArcDiagramVertex
end

function upper_vertex(::Type{ArcDiagramVertex}, i::Int)
    return (:upper, i)::ArcDiagramVertex
end

function is_upper_vertex(p::ArcDiagramVertex)
    return p[1] == :upper
end

function lower_vertices(a::ArcDiagram)
    return ArcDiagramVertex[lower_vertex(ArcDiagramVertex, i) for i in 1:n_lower_vertices(a)]
end

function lower_vertex(a::ArcDiagram, i::Int)
    @req 1 <= i <= n_lower_vertices(a) "Invalid index"
    return (:lower, i)::ArcDiagramVertex
end

function lower_vertex(::Type{ArcDiagramVertex}, i::Int)
    return (:lower, i)::ArcDiagramVertex
end

function is_lower_vertex(p::ArcDiagramVertex)
    return p[1] == :lower
end

function vertex_index(p::ArcDiagramVertex)
    return p[2]
end


################################################################################
#
# Undirected arc diagrams
#
################################################################################

struct ArcDiagramUndirected <: ArcDiagram
    n_upper_verts::Int
    n_lower_verts::Int
    upper_neighbors::Vector{ArcDiagramVertex}
    lower_neighbors::Vector{ArcDiagramVertex}

    function ArcDiagramUndirected(
        n_upper_verts::Int,
        n_lower_verts::Int,
        upper_neighbors::Vector{ArcDiagramVertex},
        lower_neighbors::Vector{ArcDiagramVertex};
        check::Bool=true,
    )
        if check
            @req length(upper_neighbors) == n_upper_verts "Upper vertices' neighbors list has wrong length."
            @req length(lower_neighbors) == n_lower_verts "Lower vertices' neighbors list has wrong length."
            for i in 1:n_upper_verts
                neigh = upper_neighbors[i]
                if is_upper_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_upper_verts "Out of bounds adjacency."
                    @req vertex_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = upper_neighbors[vertex_index(neigh)]
                    @req is_upper_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                elseif is_lower_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_lower_verts "Out of bounds adjacency."
                    neigh_neigh = lower_neighbors[vertex_index(neigh)]
                    @req is_upper_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                else
                    @req false "Invalid neighbor type."
                end
            end
            for i in 1:n_lower_verts
                neigh = lower_neighbors[i]
                if is_upper_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_upper_verts "Out of bounds adjacency."
                    neigh_neigh = upper_neighbors[vertex_index(neigh)]
                    @req is_lower_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                elseif is_lower_vertex(neigh)
                    @req 1 <= vertex_index(neigh) <= n_lower_verts "Out of bounds adjacency."
                    @req vertex_index(neigh) != i "No self-loops allowed."
                    neigh_neigh = lower_neighbors[vertex_index(neigh)]
                    @req is_lower_vertex(neigh_neigh) && vertex_index(neigh_neigh) == i "Adjacency not symmetric."
                else
                    @req false "Invalid neighbor type."
                end
            end
        end
        return new(n_upper_verts, n_lower_verts, upper_neighbors, lower_neighbors)
    end

    function ArcDiagramUndirected(
        n_upper_verts::Int,
        n_lower_verts::Int,
        upper_neighbor_inds::Vector{Int},
        lower_neighbor_inds::Vector{Int};
        check::Bool=true,
    )
        upper_neighbors = [
            i < 0 ? upper_vertex(ArcDiagramVertex, -i) : lower_vertex(ArcDiagramVertex, i) for i in upper_neighbor_inds
        ]
        lower_neighbors = [
            i < 0 ? upper_vertex(ArcDiagramVertex, -i) : lower_vertex(ArcDiagramVertex, i) for i in lower_neighbor_inds
        ]
        return ArcDiagramUndirected(n_upper_verts, n_lower_verts, upper_neighbors, lower_neighbors; check)
    end

    function ArcDiagramUndirected(upper_neighbor_inds::Vector{Int}, lower_neighbor_inds::Vector{Int}; check::Bool=true)
        n_upper_verts = length(upper_neighbor_inds)
        n_lower_verts = length(lower_neighbor_inds)
        return ArcDiagramUndirected(n_upper_verts, n_lower_verts, upper_neighbor_inds, lower_neighbor_inds; check)
    end

    function ArcDiagramUndirected(upper::AbstractString, lower::AbstractString)
        upper = strip(upper)
        lower = strip(lower)
        n_upper_verts = length(upper)
        n_lower_verts = length(lower)

        str = upper * lower
        symbols = unique(str)
        for s in symbols
            @static if VERSION >= v"1.7"
                @req count(s, str) == 2 "Symbol $s does not appear exactly twice."
            else
                @req count(string(s), str) == 2 "Symbol $s does not appear exactly twice."
            end
        end
        upper_neighbors = [(:none, 0) for _ in 1:n_upper_verts]
        lower_neighbors = [(:none, 0) for _ in 1:n_lower_verts]
        for s in symbols
            @static if VERSION >= v"1.7"
                i, j = findall(s, str)
            else
                i, j = findall(==(s), str)
            end
            if i <= n_upper_verts
                upper_neighbors[i] =
                    j <= n_upper_verts ? upper_vertex(ArcDiagramVertex, j) :
                    lower_vertex(ArcDiagramVertex, j - n_upper_verts)
            else
                lower_neighbors[i-n_upper_verts] =
                    j <= n_upper_verts ? upper_vertex(ArcDiagramVertex, j) :
                    lower_vertex(ArcDiagramVertex, j - n_upper_verts)
            end
            if j <= n_upper_verts
                upper_neighbors[j] =
                    i <= n_upper_verts ? upper_vertex(ArcDiagramVertex, i) :
                    lower_vertex(ArcDiagramVertex, i - n_upper_verts)
            else
                lower_neighbors[j-n_upper_verts] =
                    i <= n_upper_verts ? upper_vertex(ArcDiagramVertex, i) :
                    lower_vertex(ArcDiagramVertex, i - n_upper_verts)
            end
        end
        return new(n_upper_verts, n_lower_verts, upper_neighbors, lower_neighbors)
    end

    function ArcDiagramUndirected(s::AbstractString)
        upper, lower = split(s, ",")
        return ArcDiagramUndirected(upper, lower)
    end
end

function Base.:(==)(a1::ArcDiagramUndirected, a2::ArcDiagramUndirected)
    return (a1.n_upper_verts, a1.n_lower_verts, a1.upper_neighbors, a1.lower_neighbors) ==
           (a2.n_upper_verts, a2.n_lower_verts, a2.upper_neighbors, a2.lower_neighbors)
end

function Base.hash(a::ArcDiagramUndirected, h::UInt)
    b = 0x7e9086d4b4dae57d % UInt
    h = hash(a.n_upper_verts, h)
    h = hash(a.n_lower_verts, h)
    h = hash(a.upper_neighbors, h)
    h = hash(a.lower_neighbors, h)
    return xor(h, b)
end

function Base.show(io::IO, a::ArcDiagramUndirected)
    show(IOContext(io, (:compact => true)), MIME"text/plain"(), a)
end

function Base.show(io::IO, ::MIME"text/plain", a::ArcDiagramUndirected)
    symbols = join(vcat('A':'Z', 'a':'z', '0':'9'))
    if n_upper_vertices(a) + n_lower_vertices(a) > length(symbols)
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

function n_upper_vertices(a::ArcDiagramUndirected)
    return a.n_upper_verts
end

function n_lower_vertices(a::ArcDiagramUndirected)
    return a.n_lower_verts
end


function neighbor(a::ArcDiagramUndirected, v::ArcDiagramVertex)
    if is_upper_vertex(v)
        @req 1 <= vertex_index(v) <= n_upper_vertices(a) "Out of bounds access."
        return a.upper_neighbors[vertex_index(v)]
    elseif is_lower_vertex(v)
        @req 1 <= vertex_index(v) <= n_lower_vertices(a) "Out of bounds access."
        return a.lower_neighbors[vertex_index(v)]
    else
        @error "Invalid vertex."
    end
end

function neighbors(a::ArcDiagramUndirected, v::ArcDiagramVertex)
    return [neighbor(a, v)]
end

function _neighbor_of_upper_vertex(a::ArcDiagramUndirected, i::Int)
    return a.upper_neighbors[i]
end

function _neighbor_of_lower_vertex(a::ArcDiagramUndirected, i::Int)
    return a.lower_neighbors[i]
end

function is_crossing_free(a::ArcDiagramUndirected; part=:everything::Symbol)
    if part == :everything
        for v1 in vertices(a), v2 in vertices(a)
            if v1 == v2
                continue
                # symmetry of v1 and v2
            elseif is_upper_vertex(v1) == is_upper_vertex(v2) && vertex_index(v1) > vertex_index(v2)
                continue
            elseif is_lower_vertex(v1) && is_upper_vertex(v2)
                # symmetry of v1 and v2
                continue
            end
            nv1 = neighbor(a, v1)
            nv2 = neighbor(a, v2)
            if is_upper_vertex(v1) == is_upper_vertex(nv1) && vertex_index(v1) > vertex_index(nv1)
                # symmetry of v1 and nv1
                continue
            elseif is_lower_vertex(v1) && is_upper_vertex(nv1)
                # symmetry of v1 and nv1
                continue
            elseif is_upper_vertex(v2) == is_upper_vertex(nv2) && vertex_index(v2) > vertex_index(nv2)
                # symmetry of v2 and nv2
                continue
            elseif is_lower_vertex(v2) && is_upper_vertex(nv2)
                # symmetry of v2 and nv2
                continue
            elseif is_upper_vertex(v1) && is_upper_vertex(v2)
                if is_upper_vertex(nv1)
                    if is_upper_vertex(nv2)
                        if vertex_index(v2) < vertex_index(nv1) && vertex_index(nv1) < vertex_index(nv2)
                            return false
                        end
                    else
                        if vertex_index(v2) < vertex_index(nv1)
                            return false
                        end
                    end
                else
                    if is_upper_vertex(nv2)
                        continue
                    else
                        if vertex_index(nv2) < vertex_index(nv1)
                            return false
                        end
                    end
                end
            elseif is_upper_vertex(v1) && is_lower_vertex(v2)
                if is_upper_vertex(nv1)
                    continue
                else
                    if vertex_index(v2) < vertex_index(nv1) && vertex_index(nv1) < vertex_index(nv2)
                        return false
                    end
                end
            else # is_lower_vertex(v1) && is_lower_vertex(v2)
                if vertex_index(v2) < vertex_index(nv1) && vertex_index(nv1) < vertex_index(nv2)
                    return false
                end
            end
        end
        return true
    elseif part == :upper
        for v1 in upper_vertices(a), v2 in upper_vertices(a)
            # symmetry of v1 and v2
            if v1 >= v2
                continue
            end
            nv1 = neighbor(a, v1)
            nv2 = neighbor(a, v2)
            if is_lower_vertex(nv1) || is_lower_vertex(nv2)
                # not all are upper
                continue
            elseif vertex_index(v1) > vertex_index(nv1)
                # symmetry of v1 and nv1
                continue
            elseif vertex_index(v2) > vertex_index(nv2)
                # symmetry of v2 and nv2
                continue
            elseif vertex_index(v2) < vertex_index(nv1) && vertex_index(nv1) < vertex_index(nv2)
                return false
            end
        end
        return true
    elseif part == :lower
        for v1 in lower_vertices(a), v2 in lower_vertices(a)
            if v1 >= v2
                # symmetry of v1 and v2
                continue
            end
            nv1 = neighbor(a, v1)
            nv2 = neighbor(a, v2)
            if is_upper_vertex(nv1) || is_upper_vertex(nv2)
                # not all are lower
                continue
            elseif vertex_index(v1) > vertex_index(nv1)
                # symmetry of v1 and nv1
                continue
            elseif vertex_index(v2) > vertex_index(nv2)
                # symmetry of v2 and nv2
                continue
            elseif vertex_index(v2) < vertex_index(nv1) && vertex_index(nv1) < vertex_index(nv2)
                return false
            end
        end
        return true
    else
        error("Unknown part")
    end
end


struct ArcDiagramUndirectedIterator
    iter
    len::Int
end

function Base.iterate(i::ArcDiagramUndirectedIterator)
    return iterate(i.iter)
end

function Base.iterate(i::ArcDiagramUndirectedIterator, s)
    return iterate(i.iter, s)
end

Base.length(i::ArcDiagramUndirectedIterator) = i.len

Base.eltype(::Type{ArcDiagramUndirectedIterator}) = ArcDiagramUndirected


function all_arc_diagrams(
    n_upper_verts::Int,
    n_lower_verts::Int;
    indep_sets::AbstractVector{<:AbstractVector{Int}}=Vector{Int}[],
    check::Bool=true,
)
    if check
        for is in indep_sets
            @req all(i -> i < 0 ? -i <= n_upper_verts : i <= n_lower_verts, is) "Out of bounds independent sets"
        end
    end
    iter, len = iter_possible_adjacencies(
        n_upper_verts,
        n_lower_verts,
        [0 for _ in 1:n_upper_verts],
        [0 for _ in 1:n_lower_verts],
        indep_sets,
    )
    return ArcDiagramUndirectedIterator(iter, len)
end

function iter_possible_adjacencies(
    n_upper_verts::Int,
    n_lower_verts::Int,
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
            iter_possible_adjacencies(n_upper_verts, n_lower_verts, partial_upper2, partial_lower2, indep_sets)
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
                iter_possible_adjacencies(n_upper_verts, n_lower_verts, partial_upper, partial_lower2, indep_sets)
            end
            return Iterators.flatten(Iterators.map(c -> c[1], choices)), sum(c -> c[2], choices; init=0)
        else
            return [ArcDiagramUndirected(n_upper_verts, n_lower_verts, partial_upper, partial_lower; check=false)], 1
        end
    end
end
