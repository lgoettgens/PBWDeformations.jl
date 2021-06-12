struct SmashProductSymmDeformLie
    sp :: SmashProductLie
end

function Base.show(io::IO, spd::SmashProductSymmDeformLie)
    println(io, "Symmetric deformation of:")
    show(spd.sp)
end


function smashProductSymmDeformLie(dynkin::Char, n::Int64, lambda::Vector{Int64}) :: AlgebraWithCommutators{SmashProductSymmDeformLie}
    @assert n == length(lambda)
    sanitizeLieInput(dynkin, n)

    sp = smashProductLie(dynkin, n, lambda)
    basis = sp.basis
    commTable = sp.commTable
    extraData = SmashProductSymmDeformLie(sp.extraData)

    for i in 1:sp.extraData.nV, j in 1:i-1
        commTable[(mod(i), mod(j))] = []
    end

    AlgebraWithCommutators(basis, commTable, extraData)
end
