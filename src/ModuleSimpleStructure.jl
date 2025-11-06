function isomorphic_module_with_simple_structure(V::T) where {T <: LieAlgebraModule}
    if _is_standard_module(V)
        return V, identity_map(V)
    elseif ((fl, B) = _is_dual(V); fl)
        return _isomorphic_module__is_dual(V, B)
    elseif ((fl, Bs) = _is_direct_sum(V); fl)
        return _isomorphic_module__is_direct_sum(V, Bs)
    elseif ((fl, Bs) = _is_tensor_product(V); fl)
        return _isomorphic_module__is_tensor_product(V, Bs)
    elseif ((fl, B, k) = _is_exterior_power(V); fl)
        return _isomorphic_module__is_exterior_power(V, B, k)
    elseif ((fl, B, k) = _is_symmetric_power(V); fl)
        return _isomorphic_module__is_symmetric_power(V, B, k)
    elseif ((fl, B, k) = _is_tensor_power(V); fl)
        return _isomorphic_module__is_tensor_power(V, B, k)
    end
    error("not implemented for this type of module")
end

function _isomorphic_module__is_dual(V::T, B::T) where {T <: LieAlgebraModule}
    if _is_standard_module(B)
        return V, identity_map(V)
    end
    if ((fl, U) = _is_dual(B); fl)
        V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check=false)
    elseif ((fl, Cs) = _is_direct_sum(B); fl)
        U = direct_sum(dual.(Cs)...)
        V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check=false)
    elseif ((fl, Cs) = _is_tensor_product(B); fl)
        U = tensor_product(dual.(Cs)...)
        V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check=false)
    elseif ((fl, C, k) = _is_exterior_power(B); fl)
        U = exterior_power_obj(dual(C), k)
        V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check=false)
    elseif ((fl, C, k) = _is_symmetric_power(B); fl)
        inv_pure = inv(get_attribute(B, :mult_pure_function))
        U = symmetric_power_obj(dual(C), k)
        mat = zero_matrix(coefficient_ring(V), dim(V), dim(V))
        for i in 1:dim(B)
            pure_factors = inv_pure(basis(B, i))
            mat[i, i] = div(factorial(k), prod(factorial(count(==(xj), pure_factors)) for xj in unique(pure_factors)))
        end
        V_to_U = hom(V, U, mat; check=false)
    elseif ((fl, C, k) = _is_tensor_power(B); fl)
        U = tensor_power_obj(dual(C), k)
        V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check=false)
    end
    W, U_to_W = isomorphic_module_with_simple_structure(U)
    V_to_W = compose(V_to_U, U_to_W)
    return W, V_to_W
end

function _isomorphic_module__is_direct_sum(V::T, Bs::Vector{T}) where {T <: LieAlgebraModule}
    Cs_with_hom = map(isomorphic_module_with_simple_structure, Bs)
    Cs = first.(Cs_with_hom)::Vector{T}
    if length(Bs) == length(Cs) && all(splat(===), zip(Bs, Cs))
        Csum = V
        V_to_Csum = id_hom(V)
    else
        Csum = direct_sum(Cs...)
        V_to_Csum = hom_direct_sum(V, Csum, last.(Cs_with_hom)::Vector{LieAlgebraModuleHom{T, T}})
    end
    Ds = T[]
    for (C, _) in Cs_with_hom
        if ((fl, C_summands) = _is_direct_sum(C); fl)
            push!(Ds, C_summands...)
        else
            push!(Ds, C)
        end
    end
    Ds_filtered = filter(D -> dim(D) > 0, Ds)
    Ds = length(Ds_filtered) > 0 ? Ds_filtered : [Ds[1]]
    if length(Ds) == 1
        W = Ds[1]
    elseif length(Bs) == length(Ds) && all(splat(===), zip(Bs, Ds))
        W = V
    else
        W = direct_sum(Ds...)
    end
    Csum_to_W = hom(Csum, W, identity_matrix(coefficient_ring(V), dim(Csum)); check=false)
    V_to_W = compose(V_to_Csum, Csum_to_W)
    return W, V_to_W
end

function _isomorphic_module__is_tensor_product(V::T, Bs::Vector{T}) where {T <: LieAlgebraModule}
    Cs_with_hom = map(isomorphic_module_with_simple_structure, Bs)
    Cs = first.(Cs_with_hom)::Vector{T}
    if length(Bs) == length(Cs) && all(splat(===), zip(Bs, Cs))
        Cprod = V
        V_to_Cprod = id_hom(V)
    else
        Cprod = tensor_product(first.(Cs_with_hom)::Vector{T}...)
        V_to_Cprod = hom_tensor(V, Cprod, last.(Cs_with_hom)::Vector{LieAlgebraModuleHom{T, T}})
    end
    Ds = T[]
    for (C, _) in Cs_with_hom
        if ((fl, C_factors) = _is_tensor_product(C); fl)
            push!(Ds, C_factors...)
        else
            push!(Ds, C)
        end
    end
    Ds_filtered = filter(D -> dim(D) != 1 || any(!iszero, D.transformation_matrices), Ds)
    Ds = length(Ds_filtered) > 0 ? Ds_filtered : [Ds[1]]
    if length(Ds) == 1
        U = only(Ds)
    elseif length(Bs) == length(Ds) && all(splat(===), zip(Bs, Ds))
        U = V
    else
        U = tensor_product(Ds...)
    end
    Cprod_to_U = hom(Cprod, U, identity_matrix(coefficient_ring(V), dim(Cprod)); check=false)
    V_to_U = compose(V_to_Cprod, Cprod_to_U)

    if length(Ds) == 1
        return U, V_to_U
    end
    if all(D -> !_is_direct_sum(D)[1], Ds)
        W = U
        U_to_W = identity_map(U)
    else
        Es_with_summands =
            [((fl, D_summands) = _is_direct_sum(D); fl) ? (D, D_summands) : (direct_sum(D), [D]) for D in Ds]
        Fs = T[]
        inv_pure = inv(get_attribute(U, :tensor_pure_function))::MapFromFunc
        mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
        dim_accum = 0
        for summ_comb in Iterators.map(
            reverse,
            ProductIterator(reverse([1:length(E_summands) for (_, E_summands) in Es_with_summands])),
        )
            F = tensor_product([E_summands[i] for ((_, E_summands), i) in zip(Es_with_summands, summ_comb)]...)
            for (i, bi) in enumerate(basis(U))
                pure_factors = inv_pure(bi)::Tuple{Vararg{elem_type(T)}}
                dsmap = [
                    begin
                        local j, pr_f
                        Ei = Es_with_summands[i][1]
                        projs = canonical_projections(Ei)
                        if parent(f) !== Ei
                            f = Ei([f])
                        end
                        for outer j in 1:length(projs)
                            pr_f = projs[j](f)::elem_type(T)
                            if !iszero(pr_f)
                                break
                            end
                        end
                        (j, pr_f)::Tuple{Int, elem_type(T)}
                    end for (i, f::elem_type(T)) in enumerate(pure_factors)
                ]::Vector{Tuple{Int, elem_type(T)}}
                if [dsmap[l][1] for l in 1:length(Es_with_summands)] == summ_comb
                    img = F([dsmap[j][2] for j in 1:length(Es_with_summands)])
                    mat[i, dim_accum+1:dim_accum+dim(F)] = Oscar.LieAlgebras._matrix(img)
                end
            end
            dim_accum += dim(F)
            push!(Fs, F)
        end
        W = direct_sum(Fs...)
        U_to_W = hom(U, W, mat; check=false)
    end
    V_to_W = compose(V_to_U, U_to_W)
    return W, V_to_W
end

function _isomorphic_module__is_exterior_power(V::T, B::T, k::Int) where {T <: LieAlgebraModule}
    C, B_to_C = isomorphic_module_with_simple_structure(B)
    if k == 1
        U = C
        V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)); check=false)
        V_to_U = compose(V_to_B, B_to_C)
        return U, V_to_U
    end
    if B === C
        U = V
        V_to_U = id_hom(U)
    else
        U = exterior_power_obj(C, k)
        V_to_U = hom(V, U, B_to_C)
    end
    if ((fl, Ds) = _is_direct_sum(C); fl)
        m = length(Ds)
        let Ds = Ds
            Es = T[]
            inv_pure = inv(get_attribute(U, :wedge_pure_function))::MapFromFunc
            projs = canonical_projections(C)
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in multicombinations(m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? exterior_power_obj(Ds[i], lambda[i]) : nothing for i in 1:m]
                factors_cleaned = filter(!isnothing, factors)
                E = tensor_product(factors_cleaned...)
                for (i, bi) in enumerate(basis(U))
                    pure_factors = inv_pure(bi)
                    dsmap = [
                        begin
                            local j, pr_f
                            for outer j in 1:length(projs)
                                pr_f = projs[j](f)
                                if !iszero(pr_f)
                                    break
                                end
                            end
                            (j, pr_f)::Tuple{Int, elem_type(T)}
                        end for f in pure_factors
                    ]::Vector{Tuple{Int, elem_type(T)}}
                    if [dsmap[l][1] for l in 1:k] == summ_comb
                        img = E([
                            factors[j]([dsmap[l][2] for l in 1:k if dsmap[l][1] == j]) for j in 1:m if lambda[j] != 0
                        ])::elem_type(E)
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat; check=false)
            W, F_to_W = isomorphic_module_with_simple_structure(F)
            U_to_W = compose(U_to_F, F_to_W)
        end
    else
        W = U
        U_to_W = identity_map(U)
    end
    V_to_W = compose(V_to_U, U_to_W)
    return W, V_to_W
end

function _isomorphic_module__is_symmetric_power(V::T, B::T, k::Int) where {T <: LieAlgebraModule}
    C, B_to_C = isomorphic_module_with_simple_structure(B)
    if k == 1
        U = C
        V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)); check=false)
        V_to_U = compose(V_to_B, B_to_C)
        return U, V_to_U
    end
    if B === C
        U = V
        V_to_U = id_hom(U)
    else
        U = symmetric_power_obj(C, k)
        V_to_U = hom(V, U, B_to_C)
    end
    if ((fl, Ds) = _is_direct_sum(C); fl)
        m = length(Ds)
        let Ds = Ds
            Es = T[]
            inv_pure = inv(get_attribute(U, :mult_pure_function))::MapFromFunc
            projs = canonical_projections(C)
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in multicombinations(m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? symmetric_power_obj(Ds[i], lambda[i]) : nothing for i in 1:m]
                factors_cleaned = filter(!isnothing, factors)
                E = tensor_product(factors_cleaned...)
                for (i, bi) in enumerate(basis(U))
                    pure_factors = inv_pure(bi)
                    dsmap = [
                        begin
                            local j, pr_f
                            for outer j in 1:length(projs)
                                pr_f = projs[j](f)
                                if !iszero(pr_f)
                                    break
                                end
                            end
                            (j, pr_f)::Tuple{Int, elem_type(T)}
                        end for f in pure_factors
                    ]::Vector{Tuple{Int, elem_type(T)}}
                    if [dsmap[l][1] for l in 1:k] == summ_comb
                        img = E([
                            factors[j]([dsmap[l][2] for l in 1:k if dsmap[l][1] == j]) for j in 1:m if lambda[j] != 0
                        ])::elem_type(E)
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat; check=false)
            W, F_to_W = isomorphic_module_with_simple_structure(F)
            U_to_W = compose(U_to_F, F_to_W)
        end
    else
        W = U
        U_to_W = identity_map(U)
    end
    V_to_W = compose(V_to_U, U_to_W)
    return W, V_to_W
end

function _isomorphic_module__is_tensor_power(V::T, B::T, k::Int) where {T <: LieAlgebraModule}
    C, B_to_C = isomorphic_module_with_simple_structure(B)
    if k == 1
        U = C
        V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)); check=false)
        V_to_U = compose(V_to_B, B_to_C)
        return U, V_to_U
    end
    if B === C
        U = V
        V_to_U = id_hom(U)
    else
        U = tensor_power_obj(C, k)
        V_to_U = hom(V, U, B_to_C)
    end
    if ((fl, Ds) = _is_direct_sum(C); fl)
        m = length(Ds)
        let Ds = Ds
            Es = T[]
            inv_pure = inv(get_attribute(U, :tensor_pure_function))::MapFromFunc
            projs = canonical_projections(C)
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in ProductIterator(1:m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? tensor_power_obj(Ds[i], lambda[i]) : nothing for i in 1:m]
                factors_cleaned = filter(!isnothing, factors)
                E = tensor_product(factors_cleaned...)
                for (i, bi) in enumerate(basis(U))
                    pure_factors = inv_pure(bi)
                    dsmap = [
                        begin
                            local j, pr_f
                            for outer j in 1:length(projs)
                                pr_f = projs[j](f)
                                if !iszero(pr_f)
                                    break
                                end
                            end
                            (j, pr_f)::Tuple{Int, elem_type(T)}
                        end for f in pure_factors
                    ]::Vector{Tuple{Int, elem_type(T)}}
                    if [dsmap[l][1] for l in 1:k] == summ_comb
                        img = E([
                            factors[j]([dsmap[l][2] for l in 1:k if dsmap[l][1] == j]) for j in 1:m if lambda[j] != 0
                        ])::elem_type(E)
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat; check=false)
            W, F_to_W = isomorphic_module_with_simple_structure(F)
            U_to_W = compose(U_to_F, F_to_W)
        end
    else
        W = U
        U_to_W = identity_map(U)
    end
    V_to_W = compose(V_to_U, U_to_W)
    return W, V_to_W
end
