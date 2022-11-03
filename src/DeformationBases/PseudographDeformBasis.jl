function pseudographs_with_partitions__so_extpowers_stdmod(reg::Int, sumtotal::Int)
    iter = (
        begin
            (pg, Partition(copy(part)))
        end
        for sumpg in 0:sumtotal
        for pg in all_pseudographs(reg, sumpg)
        for part in AllParts(sumtotal - sumpg)
        if all(iseven, part)
        && all(isodd, pg.loops1)
        && all(isodd, pg.loops2)
        && (pg.loops1 != pg.loops2 || isodd(sum(pg.edges)))
        && pg.loops1 <= pg.loops2
    )
    return collect(iter)
end
