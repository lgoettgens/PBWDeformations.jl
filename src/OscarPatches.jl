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

if !isdefined(Oscar, :free_associative_algebra_type)
    free_associative_algebra_type(::Type{T}) where T<:RingElement = Generic.FreeAssAlgebra{T}

    free_associative_algebra_type(::Type{S}) where S<:Ring = free_associative_algebra_type(elem_type(S))
    free_associative_algebra_type(x) = free_associative_algebra_type(typeof(x)) # to stop this method from eternally recursing on itself, we better add ...
    free_associative_algebra_type(::Type{T}) where T = throw(ArgumentError("Type `$T` must be subtype of `RingElement`."))
end

#=
if !hasmethod(Oscar._is_homogeneous, Tuple{FreeAssAlgElem})
    function Oscar._is_homogeneous(f::FreeAssAlgElem)
        length(f) <= 1 && return true
        leadexpv, tailexpvs = Iterators.peel(AbstractAlgebra.exponent_words(f))
        d = length(leadexpv)
        for tailexpv in tailexpvs
        if d!=length(tailexpv)
            return false
        end
        end
        return true
    end
end
=#
