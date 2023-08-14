function isomorphic_module_with_simple_structure(V::LieAlgebraModule)
    if is_standard_module(V)
        return V, identity_map(V)
    elseif is_dual(V)
        B = base_module(V)
        if is_standard_module(B)
            return V, identity_map(V)
        end
        if is_dual(B)
            U = base_module(B)
        elseif is_direct_sum(B)
            U = direct_sum(dual.(base_modules(B))...)
        elseif is_tensor_product(B)
            U = tensor_product(dual.(base_modules(B))...)
        elseif is_exterior_power(B)
            C = base_module(B)
            k = get_attribute(B, :power)
            U = exterior_power(dual(base_module(B)), k)
        elseif is_symmetric_power(B)
            C = base_module(B)
            k = get_attribute(B, :power)
            U = symmetric_power(dual(base_module(B)), k)
        elseif is_tensor_power(B)
            C = base_module(B)
            k = get_attribute(B, :power)
            U = tensor_power(dual(base_module(B)), k)
        end
        V_to_U = hom(V, U, identity_matrix(coefficient_ring(V), dim(V)))
        W, U_to_W = isomorphic_module_with_simple_structure(U)
        V_to_W = compose(V_to_U, U_to_W)
        return W, V_to_W
    elseif is_direct_sum(V)
        Bs = base_modules(V)
        Cs_with_hom = [isomorphic_module_with_simple_structure(B) for B in Bs]
        Ds = []
        for (C, _) in Cs_with_hom
            if is_direct_sum(C)
                push!(Ds, base_modules(C)...)
            else
                push!(Ds, C)
            end
        end
        if length(Ds) == 1
            W = Ds[1]
            B_to_W = Cs_with_hom[1][2]
            h = hom(V, W, matrix(B_to_W))
            return W, h
        else
            W = direct_sum(Ds...)
            h = hom(V, W, diagonal_matrix([matrix(B_to_C) for (_, B_to_C) in Cs_with_hom]))
            return W, h
        end
    elseif is_tensor_product(V)
        Bs = base_modules(V)
        Cs_with_hom = [isomorphic_module_with_simple_structure(B) for B in Bs]
        Ds = []
        for (C, _) in Cs_with_hom
            if is_tensor_product(C)
                push!(Ds, base_modules(C)...)
            else
                push!(Ds, C)
            end
        end
        if length(Ds) == 1
            W = Ds[1]
            B_to_W = Cs_with_hom[1][2]
            h = hom(V, W, matrix(B_to_W))
            return W, h
        end
        U = tensor_product(Ds...)
        V_to_U = hom(V, U, reduce(kronecker_product, [matrix(B_to_C) for (_, B_to_C) in Cs_with_hom]))
        if all(!is_direct_sum, Ds)
            W = U
            U_to_W = identity_map(U)
        else
            Es = [is_direct_sum(D) ? D : direct_sum(D) for D in Ds]
            Fs = []
            direct_summand_mappings = [_direct_summand_mapping(E) for E in Es]
            ind_map = get_attribute(U, :ind_map)
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in reverse.(ProductIterator(reverse([1:length(base_modules(E)) for E in Es])))
                F = tensor_product([base_modules(E)[i] for (E, i) in zip(Es, summ_comb)]...)
                println(base_modules(F))
                for i in 1:dim(U)
                    inds = ind_map[i]
                    if [direct_summand_mappings[j][ind][1] for (j, ind) in enumerate(inds)] == summ_comb
                        dsmap = [direct_summand_mappings[j][ind] for (j, ind) in enumerate(inds)]
                        img = F([basis(base_modules(E)[summ_comb[j]], dsmap[j][2]) for (j, E) in enumerate(Es)])
                        mat[i, dim_accum+1:dim_accum+dim(F)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(F)
                push!(Fs, F)
            end
            W = direct_sum(Fs...)
            U_to_W = hom(U, W, mat)
        end
        V_to_W = compose(V_to_U, U_to_W)
        return W, V_to_W
    elseif is_exterior_power(V)
        B = base_module(V)
        k = get_attribute(V, :power)
        C, B_to_C = isomorphic_module_with_simple_structure(B)
        if k == 1
            U = C
            V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)))
            V_to_U = compose(V_to_B, B_to_C)
            return U, V_to_U
        end
        U = exterior_power(C, k)
        V_to_U = hom(V, U, [U(B_to_C.(_basis_repres(V, i))) for i in 1:dim(V)])
        if is_direct_sum(C)
            direct_summand_mapping = _direct_summand_mapping(C)
            ind_map = get_attribute(U, :ind_map)
            Ds = base_modules(C)
            m = length(Ds)
            Es = []
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in multicombinations(m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? exterior_power(Ds[i], lambda[i]) : nothing for i in 1:m]
                factors_cleaned = filter(!isnothing, factors)
                E = tensor_product(factors_cleaned...)
                for i in 1:dim(U)
                    inds = ind_map[i]
                    if [direct_summand_mapping[ind][1] for ind in inds] == summ_comb
                        dsmap = [direct_summand_mapping[ind] for ind in inds]
                        img = E([
                            factors[j]([basis(Ds[j], dsmap[k][2]) for k in 1:m if dsmap[k][1] == j]) for
                            j in 1:m if lambda[j] != 0
                        ])
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                if length(factors_cleaned) == 1
                    E = base_modules(E)[1]
                end
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat)
            W, F_to_W = isomorphic_module_with_simple_structure(F)
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
        C, B_to_C = isomorphic_module_with_simple_structure(B)
        if k == 1
            U = C
            V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)))
            V_to_U = compose(V_to_B, B_to_C)
            return U, V_to_U
        end
        U = symmetric_power(C, k)
        V_to_U = hom(V, U, [U(B_to_C.(_basis_repres(V, i))) for i in 1:dim(V)])
        if is_direct_sum(C)
            direct_summand_mapping = _direct_summand_mapping(C)
            ind_map = get_attribute(U, :ind_map)
            Ds = base_modules(C)
            m = length(Ds)
            Es = []
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in multicombinations(m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? symmetric_power(Ds[i], lambda[i]) : nothing for i in 1:m]
                factors_cleaned = filter(!isnothing, factors)
                E = tensor_product(factors_cleaned...)
                for i in 1:dim(U)
                    inds = ind_map[i]
                    if [direct_summand_mapping[ind][1] for ind in inds] == summ_comb
                        dsmap = [direct_summand_mapping[ind] for ind in inds]
                        img = E([
                            factors[j]([basis(Ds[j], dsmap[k][2]) for k in 1:m if dsmap[k][1] == j]) for
                            j in 1:m if lambda[j] != 0
                        ])
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                if length(factors_cleaned) == 1
                    E = base_modules(E)[1]
                end
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat)
            W, F_to_W = isomorphic_module_with_simple_structure(F)
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
        C, B_to_C = isomorphic_module_with_simple_structure(B)
        if k == 1
            U = C
            V_to_B = hom(V, B, identity_matrix(coefficient_ring(V), dim(V)))
            V_to_U = compose(V_to_B, B_to_C)
            return U, V_to_U
        end
        U = tensor_power(C, k)
        V_to_U = hom(V, U, [U(B_to_C.(_basis_repres(V, i))) for i in 1:dim(V)])
        if is_direct_sum(C)
            direct_summand_mapping = _direct_summand_mapping(C)
            ind_map = get_attribute(U, :ind_map)
            Ds = base_modules(C)
            m = length(Ds)
            Es = []
            mat = zero_matrix(coefficient_ring(U), dim(U), dim(U))
            dim_accum = 0
            for summ_comb in AbstractAlgebra.ProductIterator(1:m, k)
                lambda = [count(==(i), summ_comb) for i in 1:m]
                factors = [lambda[i] != 0 ? tensor_power(Ds[i], lambda[i]) : nothing for i in 1:m]
                factors_cleaned = filter(!isnothing, factors)
                E = tensor_product(factors_cleaned...)
                for i in 1:dim(U)
                    inds = ind_map[i]
                    if [direct_summand_mapping[ind][1] for ind in inds] == summ_comb
                        dsmap = [direct_summand_mapping[ind] for ind in inds]
                        img = E([
                            factors[j]([basis(Ds[j], dsmap[k][2]) for k in 1:m if dsmap[k][1] == j]) for
                            j in 1:m if lambda[j] != 0
                        ])
                        mat[i, dim_accum+1:dim_accum+dim(E)] = Oscar.LieAlgebras._matrix(img)
                    end
                end
                dim_accum += dim(E)
                if length(factors_cleaned) == 1
                    E = base_modules(E)[1]
                end
                push!(Es, E)
            end
            F = direct_sum(Es...)
            @assert dim(U) == dim_accum == dim(F)
            U_to_F = hom(U, F, mat)
            W, F_to_W = isomorphic_module_with_simple_structure(F)
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
    B = base_module(V)
    ind_map = get_attribute(V, :ind_map)
    js = ind_map[i]
    return map(j -> basis(B, j), js)
end

function _direct_summand_mapping(V::LieAlgebraModule)
    @req is_direct_sum(V) "Not a direct sum"
    map = [(0, 0) for _ in 1:dim(V)]
    Bs = base_modules(V)
    i = 1
    for (j, Bj) in enumerate(Bs)
        for k in 1:dim(Bj)
            map[i+k-1] = (j, k)
        end
        i += dim(Bj)
    end
    return map
end
