struct SmashProductDeformLie
    sp :: SmashProductLie
    symmetric :: Bool
    kappa :: Matrix{AlgebraElement}
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
    @assert size(kappa) == (nV, nV) "size of kappa matches module dimension"

    # basis of smash product consists of basis of module and basis of Hopf algebra
    hopfBasis = filter(!ismod, sp.basis)
    @assert all(e -> issubset(basisElements(e), hopfBasis), kappa) "kappa only takes values in Hopf algebra"

    for i in 1:nV, j in 1:i
        @assert sameSum(kappa[i,j], -kappa[j,i]) "kappa is skew-symmetric"
    end

    relTable = sp.relTable
    symmetric = true

    for i in 1:nV, j in 1:i-1
        if symmetric && kappa[i,j] != []
            symmetric = false
        end

        # We have the commutator relation [mod(i), mod(j)] = kappa[i,j]
        # which is equivalent to mod(i)*mod(j) = mod(j)*mod(i) + kappa[i,j]
        relTable[(mod(i), mod(j))] = [(1, [mod(j), mod(i)]); kappa[i,j]]
    end

    extraData = SmashProductDeformLie(sp.extraData, symmetric, kappa)
    return QuadraticAlgebra{SmashProductDeformLie}(sp.basis, relTable, extraData)
end


function smashProductSymmDeformLie(sp::QuadraticAlgebra{SmashProductLie}) :: QuadraticAlgebra{SmashProductDeformLie}
    relTable = sp.relTable

    for i in 1:sp.extraData.nV, j in 1:i-1
        relTable[(mod(i), mod(j))] = [(1, [mod(j), mod(i)])]
    end

    extraData = SmashProductDeformLie(sp.extraData, true, fill([], sp.extraData.nV, sp.extraData.nV))
    return QuadraticAlgebra{SmashProductDeformLie}(sp.basis, relTable, extraData)
end

function smashProductSymmDeformLie(dynkin::Char, n::Int64, lambda::Vector{Int64}) :: QuadraticAlgebra{SmashProductDeformLie}
    @assert n == length(lambda)
    sanitizeLieInput(dynkin, n)

    return smashProductSymmDeformLie(smashProductLie(dynkin, n, lambda))
end


function isPBWDeformation(d::QuadraticAlgebra{SmashProductDeformLie}) :: Bool
    # Uses Theorem 3.1 of Walton, Witherspoon: Poincare-Birkhoff-Witt deformations of smash product algebras from Hopf actions on Koszul algebras. DOI:	10.2140/ant.2014.8.1701. https://arxiv.org/abs/1308.6011
    nL = d.extraData.sp.nL
    nV = d.extraData.sp.nV
    kappa = d.extraData.kappa

    ## (a) κ is H-invariant
    exprs = [(sum([c*kappa[m[1][2],j] for (c,m) in normalForm(d, comm(h, mod(i)))], init=algebraElement(0)) # κ([h⋅v_i,v_j])
            + sum([c*kappa[i,m[1][2]] for (c,m) in normalForm(d, comm(h, mod(j)))], init=algebraElement(0)) # κ([v_i,h⋅v_j])
            - normalForm(d, comm(h, kappa[i,j])))                                                           # h⋅κ([v_i,v_j])
        for i in 1:nV for j in i+1:nV for h in lie(1:nL)]
    # m[1][2] denotes the index of the only basis element in the monomial m

    ## (b) trivial

    ## (c) 0 = κ ⊗ id - id ⊗ κ on (I ⊗ V) ∩ (V ⊗ I)
    # (I ⊗ V) ∩ (V ⊗ I) has basis v_iv_jv_k + v_jv_kv_i + v_kv_iv_j - v_kv_jv_i - v_jv_iv_k - v_iv_kv_j for i<j<k
    append!(exprs, [(kappa[i,j]*mod(k) - mod(i)*kappa[j,k]
            + kappa[j,k]*mod(i) - mod(j)*kappa[k,i]
            + kappa[k,i]*mod(j) - mod(k)*kappa[j,i]
            - kappa[k,j]*mod(i) + mod(k)*kappa[j,i]
            - kappa[j,i]*mod(k) + mod(j)*kappa[i,k]
            - kappa[i,k]*mod(j) + mod(i)*kappa[k,j]
        ) for i in 1:nV for j in i+1:nV for k in j+1:nV])

    ## (d) trivial

    return all(expr -> normalForm(d, expr) ≐ 0, exprs)
end
