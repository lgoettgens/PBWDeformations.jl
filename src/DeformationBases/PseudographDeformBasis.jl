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
        LieType = Val(get_attribute(base_lie_algebra(sp), :type, nothing))
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
        extra_data = Dict{DeformationMap{elem_type(sp)}, Set{Tuple{PseudographLabelled{Int}, Partition{Int}}}}()

        function diag_data_iter_and_len(LieType::SO, W::LieAlgebraModule, case::Symbol, d::Int)
            @req case == :exterior_power "Not implemented for direct sums"
            fl, V, e = _is_exterior_power(W)
            @assert fl
            @assert e == 2
            fl, Vbase, e = _is_exterior_power(V)
            @req fl && _is_standard_module(Vbase) "Only works for exterior powers of the standard module."

            pg_iter = pseudographs_with_partitions__so_extpowers_stdmod(e, d)
            diag_data_iter = ((to_arcdiag(pg, part), (pg, part)) for (pg, part) in pg_iter)
            len = length(pg_iter)
            return diag_data_iter, len::Int
        end

        function should_be_used(LieType::SO, diag::ArcDiagram, data)
            true
        end

        iter, len = arc_diag_based_basis_iteration(
            LieType,
            sp,
            degs,
            extra_data,
            diag_data_iter_and_len,
            should_be_used;
            no_normalize,
        )

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
