struct SmashProductDeformLie
    sp :: SmashProductLie
    symmetric :: Bool
end

function Base.show(io::IO, spd::SmashProductDeformLie)
    if spd.symmetric
        println(io, "Symmetric deformation of:")
    else
        println(io, "Deformation of:")
    end
    show(spd.sp)
end


function smashProductDeformLie(sp::QuadraticAlgebra{SmashProductLie}, kappa::Matrix{AlgebraElement}) :: QuadraticAlgebra{SmashProductDeformLie}
    nV = sp.extraData.nV
    @assert size(kappa) == (nV, nV)
    # @assert all(map(e -> e in sp, kappa))

    relTable = sp.relTable
    symmetric = true

    for i in 1:nV, j in 1:i-1
        if symmetric && kappa[i,j] != [(1, [mod(j), mod(i)])]
            symmetric = false
        end

        relTable[(mod(i), mod(j))] = kappa[i,j]
    end

    extraData = SmashProductDeformLie(sp.extraData, symmetric)
    QuadraticAlgebra{SmashProductDeformLie}(sp.basis, relTable, extraData)
end


function smashProductSymmDeformLie(sp::QuadraticAlgebra{SmashProductLie}) :: QuadraticAlgebra{SmashProductDeformLie}
    relTable = sp.relTable

    for i in 1:sp.extraData.nV, j in 1:i-1
        relTable[(mod(i), mod(j))] = [(1, [mod(j), mod(i)])]
    end

    extraData = SmashProductDeformLie(sp.extraData, true)
    QuadraticAlgebra{SmashProductDeformLie}(sp.basis, relTable, extraData)
end

function smashProductSymmDeformLie(dynkin::Char, n::Int64, lambda::Vector{Int64}) :: QuadraticAlgebra{SmashProductDeformLie}
    @assert n == length(lambda)
    sanitizeLieInput(dynkin, n)

    smashProductSymmDeformLie(smashProductLie(dynkin, n, lambda))
end
