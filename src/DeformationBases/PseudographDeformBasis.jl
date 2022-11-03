struct PseudographDeformBasis{C <: RingElement} <: DeformBasis{C}
    len::Int
    iter

    function PseudographDeformBasis{C}(
        sp::SmashProductLie{C},
        degs::AbstractVector{Int};
        no_normalize::Bool=false,
    ) where {C <: RingElement}
        dimV, e = extract_sp_info__so_extpowers_stdmod(sp)

        lens = []
        iters = []
        debug_counter = 0
        for d in degs
            pg_iter = pseudographs_with_partitions__so_extpowers_stdmod(e, d)
            len = length(pg_iter)
            iter = (
                begin
                    @debug "Basis generation deg $(d), $(debug_counter = (debug_counter % len) + 1)/$(len), $(floor(Int, 100*debug_counter / len))%"
                    diag = to_arcdiag(pg, part)
                    arcdiag_to_basiselem__so_extpowers_stdmod(diag, dimV, e, d, sp.alg(0), sp.basisL)
                end for (pg, part) in pg_iter
            )
            push!(lens, len)
            push!(iters, iter)
        end
        len = sum(lens)
        iter = Iterators.flatten(iters)
        if !no_normalize
            iter = normalize_basis(iter)
            len = length(iter)
        end
        return new{C}(len, iter)
    end
end

function Base.iterate(i::PseudographDeformBasis)
    return iterate(i.iter)
end

function Base.iterate(i::PseudographDeformBasis, s)
    return iterate(i.iter, s)
end

Base.length(basis::PseudographDeformBasis) = basis.len


function pseudographs_with_partitions__so_extpowers_stdmod(reg::Int, sumtotal::Int)
    iter = (
        begin
            (pg, Partition(copy(part)))
        end for sumpg in 0:sumtotal for pg in all_pseudographs(reg, sumpg) for
        part in AllParts(sumtotal - sumpg) if all(iseven, part) &&
        all(isodd, pg.loops1) &&
        all(isodd, pg.loops2) &&
        (pg.loops1 != pg.loops2 || isodd(sum(pg.edges))) &&
        pg.loops1 <= pg.loops2
    )
    return collect(iter)
end
