#################################################
#
# Abstract parent type
#
#################################################

abstract type LieAlgebraModule{C <: RingElement} <: FPModule{C} end

abstract type LieAlgebraModuleElem{C <: RingElement} <: FPModuleElem{C} end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

base_ring(v::LieAlgebraModuleElem{C}) where {C <: RingElement} = base_ring(parent(v))

gens(V::LieAlgebraModule{C}) where {C <: RingElement} = [gen(V, i) for i in 1:ngens(V)]

function gen(V::LieAlgebraModule{C}, i::Int) where {C <: RingElement}
    R = base_ring(V)
    return V([(j == i ? one(R) : zero(R)) for j in 1:ngens(V)])
end

@inline function Generic._matrix(v::LieAlgebraModuleElem{C}) where {C <: RingElement}
    return (v.mat)::dense_matrix_type(C)
end

function Generic.rels(_::LieAlgebraModule{C}) where {C <: RingElement}
    # there are no relations in a vector space
    return Vector{dense_matrix_type(C)}(undef, 0)
end

###############################################################################
#
#   String I/O
#
###############################################################################

function expressify(v::LieAlgebraModuleElem{C}, s=symbols(parent(v)); context=nothing) where {C <: RingElement}
    sum = Expr(:call, :+)
    for (i, c) in enumerate(_matrix(v))
        push!(sum.args, Expr(:call, :*, expressify(c, context=context), s[i]))
    end
    return sum
end

@enable_all_show_via_expressify LieAlgebraModuleElem


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::LieAlgebraModule{C})() where {C <: RingElement}
    mat = zero_matrix(base_ring(V), 1, ngens(V))
    return elem_type(V)(V, mat)
end

function (V::LieAlgebraModule{C})(v::Vector{C}) where {C <: RingElement}
    length(v) == ngens(V) || error("Length of vector does not match number of generators.")
    mat = matrix(base_ring(V), 1, length(v), v)
    return elem_type(V)(V, mat)
end

function (V::LieAlgebraModule{C})(v::MatElem{C}) where {C <: RingElement}
    ncols(v) == ngens(V) || error("Length of vector does not match number of generators")
    nrows(v) == 1 || error("Not a vector in module constructor")
    return elem_type(V)(V, v)
end

function (V::LieAlgebraModule{C})(v::LieAlgebraModuleElem{C}) where {C <: RingElement}
    V == parent(v) || error("Incompatible modules.")
    return v
end


###############################################################################
#
#   Arithmetic operations
#
###############################################################################

# Vector space operations get inherited from FPModule


###############################################################################
#
#   Module action
#
###############################################################################

function action(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) where {C <: RingElement}
    if haskey(parent(v).transformation_matrix_cache, x)
        transformation_matrix = parent(v).transformation_matrix_cache[x]
    else
        transformation_matrix = transformation_matrix_of_action(matrix_repr(x), parent(v))
        parent(v).transformation_matrix_cache[x] = transformation_matrix
    end
    return action_by_transformation_matrix(transformation_matrix, v)
end

function Base.:*(x::LieAlgebraElem{C}, v::LieAlgebraModuleElem{C}) where {C <: RingElement}
    return action(x, v)
end

function transformation_matrix_of_action(_::MatElem{C}, v::LieAlgebraModule{C}) where {C <: RingElement}
    error("Not implemented for $(typeof(v))")
end

function action_by_transformation_matrix(x::MatElem{C}, v::LieAlgebraModuleElem{C}) where {C <: RingElement}
    size(x, 1) == size(x, 2) || error("Transformation matrix must be square.")
    size(x, 1) == ngens(parent(v)) || error("Transformation matrix has wrong dimensions.")
    return parent(v)(_matrix(v) * transpose(x)) # equivalent to (x * v^T)^T, since we work with row vectors
end
