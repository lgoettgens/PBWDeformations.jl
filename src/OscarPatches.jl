function is_power_with_data(V::LieAlgebraModule)
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

# upstreamed in https://github.com/oscar-system/Oscar.jl/pull/4992
function Oscar.load_object(s::DeserializerState, ::Type{<:AbstractAlgebra.Generic.FreeAssociativeAlgebraElem},
                     parent_algebra::FreeAssociativeAlgebra)
  coeff_type = elem_type(base_ring(parent_algebra))
  elem = MPolyBuildCtx(parent_algebra)

  Oscar.load_array_node(s) do _
    loaded_coeff = Oscar.load_object(s, coeff_type, base_ring(parent_algebra), 2)
    loaded_term = parent_algebra(loaded_coeff)
    e = Oscar.load_array_node(s, 1) do _
      Oscar.load_object(s, Int)
    end
    # guarantees e is a Int[]
    e = convert(Vector{Int}, e)
    push_term!(elem, loaded_coeff, e)
  end

  return finish(elem)
end
