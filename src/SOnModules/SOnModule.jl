#################################################
#
# Abstract parent type
#
#################################################

abstract type SOnModule{C <: RingElement} end

abstract type SOnModuleElem{C <: RingElement} end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

# makes field access type stable
@inline function _matrix(v::SOnModuleElem{C}) where {C <: RingElement}
    return v.mat::dense_matrix_type(C)
end

base_ring(v::SOnModuleElem{C}) where {C <: RingElement} = base_ring(parent(v))

gens(V::SOnModule{C}) where {C <: RingElement} = [gen(V, i) for i in 1:ngens(V)]

function gen(V::SOnModule{C}, i::Int) where {C <: RingElement}
    R = base_ring(V)
    return V([(j == i ? one(R) : zero(R)) for j in 1:ngens(V)])
end

zero(V::SOnModule{C}) where {C <: RingElement} = V()

iszero(v::SOnModuleElem{C}) where {C <: RingElement} = iszero(_matrix(v))

function Base.hash(v::SOnModuleElem{C}, h::UInt) where {C <: RingElement}
    return hash(_matrix(v), hash(symbols(parent(v)), h))
end

function Base.deepcopy_internal(v::SOnModuleElem{C}, dict::IdDict) where {C <: RingElement}
    return parent(v)(deepcopy_internal(_matrix(v), dict))
end


###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(v::SOnModuleElem{C}, s=symbols(parent(v)); context=nothing) where {C <: RingElement}
    sum = Expr(:call, :+)
    for (i, c) in enumerate(_matrix(v))
        push!(sum.args, Expr(:call, :*, expressify(c, context=context), s[i]))
    end
    return sum
end

@enable_all_show_via_expressify SOnModuleElem


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::SOnModule{C})() where {C <: RingElement}
    mat = zero_matrix(base_ring(V), 1, ngens(V))
    return elem_type(V)(V, mat)
end

function (V::SOnModule{C})(v::Vector{C}) where {C <: RingElement}
    length(v) == ngens(V) || error("Length of vector does not match number of generators.")
    mat = matrix(base_ring(V), 1, length(v), v)
    return elem_type(V)(V, mat)
end

function (V::SOnModule{C})(v::MatElem{C}) where {C <: RingElement}
    ncols(v) == ngens(V) || error("Length of vector does not match number of generators")
    nrows(v) == 1 || error("Not a vector in module constructor")
    return elem_type(V)(V, v)
end

function (V::SOnModule{C})(v::SOnModuleElem{C}) where {C <: RingElement}
    V == parent(v) || error("Incompatible modules.")
    return v
end


###############################################################################
#
#   Arithmetic operations
#
###############################################################################

function Base.:-(v::SOnModuleElem{C}) where {C <: RingElement}
    return parent(v)(-_matrix(v))
end

function Base.:+(v1::SOnModuleElem{C}, v2::SOnModuleElem{C}) where {C <: RingElement}
    parent(v1) == parent(v2) || error("Incompatible modules.")
    return parent(v1)(_matrix(v1) + _matrix(v2))
end

function Base.:-(v1::SOnModuleElem{C}, v2::SOnModuleElem{C}) where {C <: RingElement}
    parent(v1) == parent(v2) || error("Incompatible modules.")
    return parent(v1)(_matrix(v1) - _matrix(v2))
end

function Base.:*(c::C, v::SOnModuleElem{C}) where {C <: RingElement}
    return parent(v)(c * _matrix(v))
end

function Base.:*(c::Union{Integer, Rational, AbstractFloat}, v::SOnModuleElem{C}) where {C <: RingElement}
    return parent(v)(c * _matrix(v))
end

Base.:*(v::SOnModuleElem{C}, c::Union{Integer, Rational, AbstractFloat}) where {C <: RingElement} = c * v

Base.:*(v::SOnModuleElem{C}, c::C) where {C <: RingElement} = c * v


###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(v1::SOnModuleElem{C}, v2::SOnModuleElem{C}) where {C <: RingElement}
    parent(v1) == parent(v2) || return false
    return _matrix(v1) == _matrix(v2)
end