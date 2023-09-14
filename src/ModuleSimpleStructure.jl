function isomorphic_module_with_simple_structure(V::LieAlgebraModule; check::Bool=true)
    if is_standard_module(V)
        return V, identity_map(V)
    elseif is_dual(V)
        B = base_module(V)
        if is_standard_module(B)
            return V, identity_map(V)
        end
        if is_dual(B)
            U = base_module(B)
            V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check)
        elseif is_direct_sum(B)
            U = direct_sum(dual.(base_modules(B))...)
            V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check)
        elseif is_tensor_product(B)
            U = tensor_product(dual.(base_modules(B))...)
            V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check)
        elseif is_exterior_power(B)
            C = base_module(B)
            k = get_attribute(B, :power)
            U = exterior_power(dual(base_module(B)), k)
            V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check)
        elseif is_symmetric_power(B)
            C = base_module(B)
            k = get_attribute(B, :power)
            inv_pure = get_attribute(B, :symmetric_pure_preimage_function)
            U = symmetric_power(dual(base_module(B)), k)
            mat = zero_matrix(coefficient_ring(V), dim(V), dim(V))
            for i in 1:dim(B)
                pure_factors = inv_pure(basis(B, i))
                mat[i, i] =
                    div(factorial(k), prod(factorial(count(==(xj), pure_factors)) for xj in unique(pure_factors)))
            end
            V_to_U = hom(V, U, mat; check)
        elseif is_tensor_power(B)
            C = base_module(B)
            k = get_attribute(B, :power)
            U = tensor_power(dual(base_module(B)), k)
            V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)); check)
        end
        W, U_to_W = isomorphic_module_with_simple_structure(U; check)
        V_to_W = compose(V_to_U, U_to_W)
        return W, V_to_W
    elseif is_direct_sum(V)
        Bs = base_modules(V)
        Cs_with_hom = [isomorphic_module_with_simple_structure(B; check) for B in Bs]
        Ds = []
        for (C, _) in Cs_with_hom
            if is_direct_sum(C)
                push!(Ds, base_modules(C)...)
            else
                push!(Ds, C)
            end
        end
        Ds_filtered = filter(D -> dim(D) > 0, Ds)
        Ds = length(Ds_filtered) > 0 ? Ds_filtered : [Ds[1]]
        if length(Ds) == 1
            W = Ds[1]
            i = findfirst(C_with_hom -> C_with_hom[1] == W, Cs_with_hom)
            B_to_W = Cs_with_hom[i][2]
            h = hom(V, W, matrix(B_to_W); check)
            return W, h
        else
            W = direct_sum(Ds...)
            h = hom(V, W, diagonal_matrix([matrix(B_to_C) for (_, B_to_C) in Cs_with_hom]); check)
            return W, h
        end
    elseif is_tensor_product(V)
        Bs = base_modules(V)
        Cs_with_hom = [isomorphic_module_with_simple_structure(B; check) for B in Bs]
        Ds = []
        for (C, _) in Cs_with_hom
            if is_tensor_product(C)
                push!(Ds, base_modules(C)...)
            else
                push!(Ds, C)
            end
        end
        Ds_filtered = filter(D -> dim(D) != 1 || any(!iszero, D.transformation_matrices), Ds)
        Ds = length(Ds_filtered) > 0 ? Ds_filtered : [Ds[1]]
        if length(Ds) == 1
            W = Ds[1]
            i = findfirst(C_with_hom -> C_with_hom[1] == W, Cs_with_hom)
            B_to_W = Cs_with_hom[i][2]
            h = hom(V, W, matrix(B_to_W); check)
            return W, h
        end
        U = tensor_product(Ds...)
        V_to_U = hom(V, U, reduce(kronecker_product, [matrix(B_to_C) for (_, B_to_C) in Cs_with_hom]); check)
        if all(!is_direct_sum, Ds)
            W = U
            U_to_W = identity_map(U)
        else
            Es = [is_direct_sum(D) ? D : direct_sum(D) for D in Ds]
            Fs = []
            inv_pure = get_attribute(U, :tensor_pure_preimage_function)
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in reverse.(ProductIterator(reverse([1:length(base_modules(E)) for E in Es])))
                F = tensor_product([base_modules(E)[i] for (E, i) in zip(Es, summ_comb)]...)
                for (i, bi) in enumerate(basis(U))
                    pure_factors = inv_pure(bi)
                    dsmap = [
                        begin
                            local j, pr_f
                            projs = canonical_projections(Es[i])
                            if parent(f) !== Es[i]
                                f = Es[i]([f])
                            end
                            for outer j in 1:length(projs)
                                pr_f = projs[j](f)
                                if !iszero(pr_f)
                                    break
                                end
                            end
                            (j, pr_f)
                        end for (i, f) in enumerate(pure_factors)
                    ]
                    if [dsmap[l][1] for l in 1:length(Es)] == summ_comb
                        img = F([dsmap[j][2] for (j, E) in enumerate(Es)])
                        mat[i, dim_accum+1:dim_accum+dim(F)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(F)
                push!(Fs, F)
            end
            W = direct_sum(Fs...)
            U_to_W = hom(U, W, mat; check)
        end
        V_to_W = compose(V_to_U, U_to_W)
        return W, V_to_W
    elseif is_exterior_power(V)
        B = base_module(V)
        k = get_attribute(V, :power)
        C, B_to_C = isomorphic_module_with_simple_structure(B; check)
        if k == 1
            U = C
            V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)); check)
            V_to_U = compose(V_to_B, B_to_C)
            return U, V_to_U
        end
        U = exterior_power(C, k)
        V_to_U = hom(V, U, elem_type(U)[U(B_to_C.(_basis_repres(V, i))) for i in 1:dim(V)]; check)
        if is_direct_sum(C)
            Ds = base_modules(C)
            m = length(Ds)
            Es = []
            inv_pure = get_attribute(U, :exterior_pure_preimage_function)
            projs = canonical_projections(C)
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in multicombinations(m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? exterior_power(Ds[i], lambda[i]) : nothing for i in 1:m]
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
                            (j, pr_f)
                        end for f in pure_factors
                    ]
                    if [dsmap[l][1] for l in 1:k] == summ_comb
                        img = E([
                            factors[j]([dsmap[l][2] for l in 1:k if dsmap[l][1] == j]) for j in 1:m if lambda[j] != 0
                        ])
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat; check)
            W, F_to_W = isomorphic_module_with_simple_structure(F; check)
            U_to_W = compose(U_to_F, F_to_W)
        else
            W = U
            U_to_W = identity_map(U)
        end
        V_to_W = compose(V_to_U, U_to_W)
        return W, V_to_W
    elseif is_symmetric_power(V)
        B = base_module(V)
        k = get_attribute(V, :power)
        C, B_to_C = isomorphic_module_with_simple_structure(B; check)
        if k == 1
            U = C
            V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)); check)
            V_to_U = compose(V_to_B, B_to_C)
            return U, V_to_U
        end
        U = symmetric_power(C, k)
        V_to_U = hom(V, U, elem_type(U)[U(B_to_C.(_basis_repres(V, i))) for i in 1:dim(V)]; check)
        if is_direct_sum(C)
            Ds = base_modules(C)
            m = length(Ds)
            Es = []
            inv_pure = get_attribute(U, :symmetric_pure_preimage_function)
            projs = canonical_projections(C)
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in multicombinations(m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? symmetric_power(Ds[i], lambda[i]) : nothing for i in 1:m]
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
                            (j, pr_f)
                        end for f in pure_factors
                    ]
                    if [dsmap[l][1] for l in 1:k] == summ_comb
                        img = E([
                            factors[j]([dsmap[l][2] for l in 1:k if dsmap[l][1] == j]) for j in 1:m if lambda[j] != 0
                        ])
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat; check)
            W, F_to_W = isomorphic_module_with_simple_structure(F; check)
            U_to_W = compose(U_to_F, F_to_W)
        else
            W = U
            U_to_W = identity_map(U)
        end
        V_to_W = compose(V_to_U, U_to_W)
        return W, V_to_W
    elseif is_tensor_power(V)
        B = base_module(V)
        k = get_attribute(V, :power)
        C, B_to_C = isomorphic_module_with_simple_structure(B; check)
        if k == 1
            U = C
            V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)); check)
            V_to_U = compose(V_to_B, B_to_C)
            return U, V_to_U
        end
        U = tensor_power(C, k)
        V_to_U = hom(V, U, elem_type(U)[U(B_to_C.(_basis_repres(V, i))) for i in 1:dim(V)]; check)
        if is_direct_sum(C)
            Ds = base_modules(C)
            m = length(Ds)
            Es = []
            inv_pure = get_attribute(U, :tensor_pure_preimage_function)
            projs = canonical_projections(C)
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in AbstractAlgebra.ProductIterator(1:m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? tensor_power(Ds[i], lambda[i]) : nothing for i in 1:m]
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
                            (j, pr_f)
                        end for f in pure_factors
                    ]
                    if [dsmap[l][1] for l in 1:k] == summ_comb
                        img = E([
                            factors[j]([dsmap[l][2] for l in 1:k if dsmap[l][1] == j]) for j in 1:m if lambda[j] != 0
                        ])
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat; check)
            W, F_to_W = isomorphic_module_with_simple_structure(F; check)
            U_to_W = compose(U_to_F, F_to_W)
        else
            W = U
            U_to_W = identity_map(U)
        end
        V_to_W = compose(V_to_U, U_to_W)
        return W, V_to_W
    end
    error("not implemented for this type of module")
end


function _basis_repres(V::LieAlgebraModule, i::Int)
    @req is_exterior_power(V) || is_symmetric_power(V) || is_tensor_power(V) "Not a power module."
    inv_pure = if is_exterior_power(V)
        get_attribute(V, :exterior_pure_preimage_function)
    elseif is_symmetric_power(V)
        get_attribute(V, :symmetric_pure_preimage_function)
    elseif is_tensor_power(V)
        get_attribute(V, :tensor_pure_preimage_function)
    end
    return inv_pure(basis(V, i))
end
