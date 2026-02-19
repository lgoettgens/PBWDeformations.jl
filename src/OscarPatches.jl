function is_power_with_data(V::LieAlgebraModuleOrLazy)
    data = _is_exterior_power(V)
    data[1] && return data
    data = _is_symmetric_power(V)
    data[1] && return data
    data = _is_tensor_power(V)
    return data
end


function exterior_power_obj(V::LieAlgebraModule, k::Int)
    return exterior_power(V, k)[1]
end

function symmetric_power_obj(V::LieAlgebraModule, k::Int)
    return symmetric_power(V, k)[1]
end

function tensor_power_obj(V::LieAlgebraModule, k::Int)
    return tensor_power(V, k)[1]
end

function is_smallest_obj_in_orbit(g::T, G::PermGroup; lt::Function=_lt) where T
    is_fixpoint = true
    for p in gens(G)
        gp = g ^ p
        is_fixpoint = is_fixpoint && gp == g
        lt(gp, g) && return false
    end
    is_fixpoint && return true
    return !any(gp -> lt(gp, g), orbit(G, ^, g))
end
