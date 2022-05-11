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
    rel_indep_sets = filter(is -> in(i, is), indep_sets)
    poss_adjs = setdiff(setdiff(findall(==(0), partial), i), rel_indep_sets...)
    choices = Iterators.map(poss_adjs) do j
        partial2 = deepcopy(partial)
        partial2[i] = j
        partial2[j] = i
        return iter_possible_adjacencies(nUpper, nLower, indep_sets, partial2)
    end
    return Iterators.flatten(Iterators.map(c -> c[1], choices)), sum(c -> c[2], choices; init=0)
end

function pbw_arc_diagrams(l::Int, d::Int)
    indep_sets = [1:l, l+1:2*l, [[2 * l + 2 * i - 1, 2 * l + 2 * i] for i in 1:d]...]
    return all_arc_diagrams(2 * l, 2 * d; indep_sets)
end

struct SoDeformArcBase{C <: RingElement} <: DeformBase{C}
    len::Int
    iter

    function SoDeformArcBase{C}(sp::SmashProductLie{C}, maxdeg::Int) where {C <: RingElement}
        dimV = div(isqrt(8 * sp.dimL + 1) + 1, 2) # inverse of gaussian sum
        l = first(l for l in 1:dimV if binomial(dimV, l) >= sp.dimV)

        d = maxdeg
        iso_wedge2V_g = Dict{Vector{Int}, Int}()
        for (i, bs) in enumerate(Combinatorics.combinations(1:dimV, 2))
            iso_wedge2V_g[bs] = i
        end
        index = Dict{Vector{Int}, Int}()
        for (i, is) in enumerate(Combinatorics.combinations(1:dimV, l))
            index[is] = i
        end

        diag_iter = pbw_arc_diagrams(l, d)
        len = length(diag_iter)
        iter = (
            begin
                kappa = fill(sp.alg(0), sp.dimV, sp.dimV)
                for is in Combinatorics.combinations(1:dimV, l), js in Combinatorics.combinations(1:dimV, l)
                    i = index[is]
                    j = index[js]
                    if i >= j
                        continue
                    end
                    zeroprod = false
                    labeled_diag = [is..., js..., [0 for _ in 1:2d]...]
                    frees = Int[]
                    for k in 1:length(labeled_diag)
                        if labeled_diag[k] != 0
                            if labeled_diag[diag.adjacency[k]] == 0
                                labeled_diag[diag.adjacency[k]] = labeled_diag[k]
                            else
                                if labeled_diag[k] != labeled_diag[diag.adjacency[k]]
                                    zeroprod = true
                                    break
                                end
                            end
                        else
                            if labeled_diag[diag.adjacency[k]] == 0
                                append!(frees, min(k, diag.adjacency[k]))
                            end
                        end
                    end
                    if zeroprod
                        continue
                    end
                    unique!(sort!(frees))
                    free_index = Dict{Int, Int}()
                    for (k, f) in enumerate(frees)
                        free_index[f] = k
                    end
                    entry = zero(sp.alg)
                    for labeling in (isempty(frees) ? [Int[]] : AbstractAlgebra.ProductIterator(1:dimV, length(frees)))
                        lower_labeled = [
                            labeled_diag[k] != 0 ? labeled_diag[k] : labeling[free_index[min(k, diag.adjacency[k])]] for k in 2*l+1:2*l+2*d
                        ]
                        zeroelem = false
                        sign_pos = true
                        basiselem = Int[]
                        for k in 1:2:length(lower_labeled)
                            if lower_labeled[k] == lower_labeled[k+1]
                                zeroelem = true
                                break
                            elseif lower_labeled[k] > lower_labeled[k+1]
                                sign_pos = !sign_pos
                                append!(basiselem, iso_wedge2V_g[[lower_labeled[k+1], lower_labeled[k]]])
                            else
                                append!(basiselem, iso_wedge2V_g[[lower_labeled[k], lower_labeled[k+1]]])
                            end
                        end
                        if zeroelem
                            continue
                        end
                        symm_basiselem =
                            1 // factorial(length(basiselem)) *
                            sum(prod(sp.baseL[ind]) for ind in Combinatorics.permutations(basiselem))
                        entry += (sign_pos ? 1 : (-1)) * normal_form(symm_basiselem)
                    end
                    kappa[i, j] += entry
                    kappa[j, i] -= entry
                end
                iszero(kappa) ? nothing : kappa
            end for diag in diag_iter
        )
        return new{C}(len, iter)
    end
end

function Base.iterate(i::SoDeformArcBase)
    return iterate(i.iter)
end

function Base.iterate(i::SoDeformArcBase, s)
    return iterate(i.iter, s)
end

Base.length(base::SoDeformArcBase) = base.len
