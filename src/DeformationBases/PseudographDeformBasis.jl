"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by a pseudograph with two vertices and
certain properties, which gets transformed to an arc diagram and then handled as
in [`ArcDiagDeformBasis`](@ref).
This process is due to [FM22](@cite).
"""
struct PseudographDeformBasis{T <: SmashProductLieElem} <: DeformBasis{T}
    len::Int
    iter
    extra_data::Dict{DeformationMap{T}, Set{Tuple{PseudographLabelled{Int}, Partition{Int}}}}
    no_normalize::Bool

    function PseudographDeformBasis(
        sp::SmashProductLie{C, LieC, LieT},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        LieType = Val(get_attribute(base_lie_algebra(sp), :type, nothing)::Union{Nothing, Symbol})
        @req LieType isa SO "Only works for so_n."
        if LieType isa SO && has_attribute(base_lie_algebra(sp), :form)
            @req isone(get_attribute(base_lie_algebra(sp), :form)) "Only works for so_n represented as skew-symmetric matrices."
        end
        return PseudographDeformBasis(LieType, sp, degs; no_normalize)
    end

    function PseudographDeformBasis(
        LieType::SO,
        sp::SmashProductLie{C, LieC, LieT},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
        fl, Vbase, e = _is_exterior_power(base_module(sp))
        @req fl && _is_standard_module(Vbase) "Only works for exterior powers of the standard module."

        extra_data = Dict{DeformationMap{elem_type(sp)}, Set{Tuple{PseudographLabelled{Int}, Partition{Int}}}}()

        lens = []
        iters = []
        debug_counter = 0
        for d in degs
            pg_iter = pseudographs_with_partitions__so_extpowers_stdmod(e, d)
            len = length(pg_iter)
            iter = (
                begin
                    @vprintln :PBWDeformations 2 "Basis generation deg $(lpad(d, maximum(ndigits, degs))), $(lpad(floor(Int, 100*(debug_counter = (debug_counter % len) + 1) / len), 3))%, $(lpad(debug_counter, ndigits(len)))/$(len)"
                    diag = to_arcdiag(pg, part)
                    basis_elem = arcdiag_to_deformationmap(LieType, diag, sp)
                    if !no_normalize
                        basis_elem = normalize(basis_elem)
                    end
                    if haskey(extra_data, basis_elem)
                        push!(extra_data[basis_elem], (pg, part))
                    else
                        extra_data[basis_elem] = Set([(pg, part)])
                    end
                    basis_elem
                end for (pg, part) in pg_iter
            )
            push!(lens, len)
            push!(iters, iter)
        end
        len = sum(lens)
        iter = Iterators.flatten(iters)
        if !no_normalize
            iter = unique(Iterators.filter(b -> !iszero(b), iter))
            collected = Vector{DeformationMap{elem_type(sp)}}(collect(iter))::Vector{DeformationMap{elem_type(sp)}}
            _, rels = is_linearly_independent_with_relations(coefficient_ring(sp), reverse(collected))
            inds = [1 + ncols(rels) - (findfirst(!iszero, vec(rels[i, :]))::Int) for i in nrows(rels):-1:1] # FIXME: findfirst -> findlast
            deleteat!(collected, inds)
            return new{elem_type(sp)}(length(collected), collected, extra_data, no_normalize)
        end
        return new{elem_type(sp)}(len, iter, extra_data, no_normalize)
    end
end

function Base.iterate(i::PseudographDeformBasis)
    return iterate(i.iter)::Union{Tuple{eltype(i), Any}, Nothing}
end

function Base.iterate(i::PseudographDeformBasis, s)
    return iterate(i.iter, s)::Union{Tuple{eltype(i), Any}, Nothing}
end

Base.length(basis::PseudographDeformBasis) = basis.len


function pseudographs_with_partitions__so_extpowers_stdmod(deg::Int, sumtotal::Int)
    iter = (
        begin
            (pg, part)
        end for sumpg in 0:sumtotal for pg in all_pseudographs(2, deg, sumpg; upto_iso=true) for
        part in partitions(sumtotal - sumpg) if all(iseven, part) &&
        all(isodd, edge_labels(pg, MSet([1, 1]))) &&
        all(isodd, edge_labels(pg, MSet([2, 2]))) &&
        (edges(pg, MSet([1, 1])) != edges(pg, MSet([1, 1])) || isodd(sum(pg, MSet([1, 2]))))
    )
    return collect(iter)
end
