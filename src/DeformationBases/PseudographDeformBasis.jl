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
        extra_data = Dict{DeformationMap{elem_type(sp)}, Set{Tuple{PseudographLabelled{Int}, Partition{Int}}}}()

        function diag_data_iter_and_len(LieType::SO, W::LieAlgebraModule, case::Symbol, d::Int)
            @req case == :exterior_power "Not implemented for direct sums"
            fl, V, e = _is_exterior_power(W)
            @assert fl
            @assert e == 2
            fl, Vbase, e = _is_exterior_power(V)
            @req fl && _is_standard_module(Vbase) "Only works for exterior powers of the standard module."

            pg_part_iter = collect(
                begin
                    (pg, part)
                end for sumpg in 0:d for pg in all_pseudographs(2, e, sumpg; upto_iso=true) for
                part in partitions(d - sumpg)

            )
            diag_data_iter = ((to_arcdiag(pg, part), (pg, part)) for (pg, part) in pg_part_iter)
            len = length(pg_part_iter) # TODO: implement without collect
            return diag_data_iter, len::Int
        end

        function should_be_used(LieType::SO, diag::ArcDiagram, data::Tuple{PseudographLabelled{Int}, Partition{Int}})
            pg, part = data

            # length n lower cycles can be reversed in n lower flips of sgn -1 => odd n has sgn -1 in stabilizer
            all(iseven, part) || return false

            # length n paths from 1 to 1 can be reversed with n-1 many lower flips of sgn -1 and one upper flip of sgn 1 => even n has sgn -1 in stabilizer
            all(isodd, edge_labels(pg, MSet([1, 1]))) || return false

            # even paths from 2 to 2 can be reversed with odd many lower flips of sgn -1 and one upper flip of sgn 1 => even n has sgn -1 in stabilizer
            all(isodd, edge_labels(pg, MSet([2, 2]))) || return false

            # if pg is symmetric, this symmetry needs one upper block flip of sgn -1 and sum(pg, MSet([1, 2])) many lower flips of sgn -1 => even sum has sgn -1 in stabilizer
            (edges(pg, MSet([1, 1])) != edges(pg, MSet([2, 2])) || isodd(sum(pg, MSet([1, 2])))) || return false

            return true
        end

        iter1, len1 = arc_diag_based_basis_iteration(
            LieType,
            sp,
            degs,
            extra_data,
            diag_data_iter_and_len,
            should_be_used;
            no_normalize,
        )

        if no_normalize
            return new{elem_type(sp)}(len1, iter1, extra_data, no_normalize)
        else
            iter2, len2 = filter_independent(coefficient_ring(sp), iter1)
            return new{elem_type(sp)}(len2, iter2, extra_data, no_normalize)
        end
    end
end

function Base.iterate(i::PseudographDeformBasis)
    return iterate(i.iter)::Union{Tuple{eltype(i), Any}, Nothing}
end

function Base.iterate(i::PseudographDeformBasis, s)
    return iterate(i.iter, s)::Union{Tuple{eltype(i), Any}, Nothing}
end

Base.length(basis::PseudographDeformBasis) = basis.len
