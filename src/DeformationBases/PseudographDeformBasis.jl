const PseudographDeformBasisDataT = Tuple{PseudographLabelled{Int}, Partition{Int}}

"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by a pseudograph with two vertices and
certain properties, which gets transformed to an arc diagram and then handled as
in [`ArcDiagDeformBasis`](@ref).
This process is due to [FM22](@cite).
"""
const PseudographDeformBasis{T} =
    ArcDiagBasedDeformBasis{PseudographDeformBasisDataT, T} where {T <: SmashProductLieElem}

function check_input(
    ::Type{PseudographDeformBasis},
    LieType,
    sp::SmashProductLie{C, LieC, LieT},
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    @req LieType isa SO "Only works for so_n."
    @req has_attribute(base_lie_algebra(sp), :form) && isone(get_attribute(base_lie_algebra(sp), :form)::dense_matrix_type(LieC)) "Only works for so_n represented as skew-symmetric matrices."
end

function data_iter_and_len(::Type{PseudographDeformBasis}, LieType::SO, W::LieAlgebraModule, case::Symbol, d::Int)
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
    len = length(pg_part_iter) # TODO: implement without collect
    return pg_part_iter, len::Int
end

function should_data_be_used(
    ::Type{PseudographDeformBasis},
    LieType::SO,
    data::PseudographDeformBasisDataT,
    ::SmashProductLie,
    ::LieAlgebraModule,
    ::Symbol,
    cache::Union{Dict{<:Any, Bool}, Nothing},
)
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

