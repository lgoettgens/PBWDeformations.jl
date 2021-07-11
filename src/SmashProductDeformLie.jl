struct SmashProductDeformLie
    sp :: SmashProductLie
    symmetric :: Bool
end

function Base.:(==)(spd1::SmashProductDeformLie, spd2::SmashProductDeformLie) :: Bool
    (spd1.sp, spd1.symmetric) == (spd2.sp, spd2.symmetric)
end

function Base.show(io::IO, spd::SmashProductDeformLie) :: Nothing
    if spd.symmetric
        println(io, "Symmetric deformation of:")
    else
        println(io, "Deformation of:")
    end
    print(io, spd.sp)
end


function smashProductDeformLie(sp::QuadraticAlgebra{SmashProductLie}, kappa::Matrix{AlgebraElement}) :: QuadraticAlgebra{SmashProductDeformLie}
    nV = sp.extraData.nV
    @assert size(kappa) == (nV, nV) "size of kappa does not match module dimension"

    # basis of smash product consists of basis of module and basis of Hopf algebra
    hopfBasis = filter(!isMod, sp.basis)
    @assert all(e -> issubset(basisElements(e), hopfBasis), kappa) "kappa does not only take values in Hopf algebra"

    relTable = sp.relTable
    symmetric = true

    for i in 1:nV, j in 1:i-1
        #@assert kappa[i,j] == -kappa[j,i] "kappa is not skew-symmetric"

        if symmetric && kappa[i,j] != []
            symmetric = false
        end

        # We have the commutator relation [mod(i), mod(j)] = kappa[i,j]
        # which is equivalent to mod(i)*mod(j) = mod(j)*mod(i) + kappa[i,j]
        relTable[(mod(i), mod(j))] = [(1, [mod(j), mod(i)]); kappa[i,j]]
    end

    extraData = SmashProductDeformLie(sp.extraData, symmetric)
    return QuadraticAlgebra{SmashProductDeformLie}(sp.basis, relTable, extraData)
end


function smashProductSymmDeformLie(sp::QuadraticAlgebra{SmashProductLie}) :: QuadraticAlgebra{SmashProductDeformLie}
    relTable = sp.relTable

    for i in 1:sp.extraData.nV, j in 1:i-1
        relTable[(mod(i), mod(j))] = [(1, [mod(j), mod(i)])]
    end

    extraData = SmashProductDeformLie(sp.extraData, true)
    return QuadraticAlgebra{SmashProductDeformLie}(sp.basis, relTable, extraData)
end

function smashProductSymmDeformLie(dynkin::Char, n::Int64, lambda::Vector{Int64}) :: QuadraticAlgebra{SmashProductDeformLie}
    @assert n == length(lambda)
    sanitizeLieInput(dynkin, n)

    return smashProductSymmDeformLie(smashProductLie(dynkin, n, lambda))
end
