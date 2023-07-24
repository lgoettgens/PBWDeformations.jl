"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by a pseudograph with two vertices and
certain properties, which gets transformed to an arc diagram and then handled as
in [`ArcDiagDeformBasis`](@ref).
This process is due to [FM22](@cite).
"""
struct PseudographDeformBasis{C <: RingElem} <: DeformBasis{C}
    len::Int
    iter
    extra_data::Dict{DeformationMap{C}, Set{Tuple{PseudographLabelled{Int}, Partition{Int}}}}
    normalize

    function PseudographDeformBasis{C}(
        sp::SmashProductLie{C},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElem}
        @req get_attribute(base_lie_algebra(sp), :type, nothing) == :special_orthogonal "Only works for so_n."
        @req is_exterior_power(base_module(sp)) && is_standard_module(base_module(base_module(sp))) "Only works for exterior powers of the standard module."

        e = get_attribute(base_module(sp), :power)

        extra_data = Dict{DeformationMap{C}, Set{Tuple{PseudographLabelled{Int}, Partition{Int}}}}()
        normalize = no_normalize ? identity : normalize_default

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
                    basis_elem = arcdiag_to_deformationmap__so(diag, sp)
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
            len = length(iter)
        end
        return new{C}(len, iter, extra_data, normalize)
    end
end

function Base.iterate(i::PseudographDeformBasis)
    return iterate(i.iter)
end

function Base.iterate(i::PseudographDeformBasis, s)
    return iterate(i.iter, s)
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
