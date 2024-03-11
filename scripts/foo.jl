function printboth(io::IO, xs...)
    print(io, xs...)
    flush(io)
    print(xs...)
end

function foo(L, V, d)
    path, fileio = mktemp()
    println(path)
    show(fileio, MIME"text/plain"(), L)
    print(fileio, "\n\n")
    show(fileio, MIME"text/plain"(), V)
    print(fileio, "\n\n")
    printboth(fileio, "deg", "\t")
    printboth(fileio, "new", "\t")
    printboth(fileio, "total", "\n")
    try
        sp = smash_product(L, V)
        bs = ArcDiagDeformBasis{QQFieldElem}[]
        for i in 0:d
            printboth(fileio, i, "\t")
            push!(bs, ArcDiagDeformBasis{QQFieldElem}(sp, i:i))
            ms = all_pbwdeformations(sp, bs[end])
            printboth(fileio, length(ms), "\t")
            ms = all_pbwdeformations(sp, UnionBasis{QQFieldElem}(bs))
            printboth(fileio, length(ms), "\n")
        end
        println("Done.")
    catch e
        if e isa InterruptException
            println("Interrupted.")
        else
            rethrow(e)
        end
    finally
        close(fileio)
    end
    println("Total results:")
    print(read(path, String))
end
