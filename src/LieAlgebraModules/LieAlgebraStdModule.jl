struct LieAlgebraStdModule{C <: RingElement} <: LieAlgebraModule{C}
    R::Ring
    n::Int
    transformation_matrix_cache::Dict{MatElem{C}, MatElem{C}}

    function LieAlgebraStdModule{C}(R::Ring, n::Int) where {C <: RingElement}
        transformation_matrix_cache = Dict{MatElem{C}, MatElem{C}}()
        return new{C}(R, n, transformation_matrix_cache)
    end
end

struct LieAlgebraStdModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraStdModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraStdModuleElem{C}}) where {C <: RingElement} = LieAlgebraStdModule{C}

elem_type(::Type{LieAlgebraStdModule{C}}) where {C <: RingElement} = LieAlgebraStdModuleElem{C}

parent(v::LieAlgebraStdModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraStdModule{C}) where {C <: RingElement} = V.R

ngens(V::LieAlgebraStdModule{C}) where {C <: RingElement} = V.n


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraStdModule{C}) where {C <: RingElement}
    print(io, "SO_$(V.n)_StdModule over ")
    print(IOContext(io, :compact => true), base_ring(V))
end

function symbols(V::LieAlgebraStdModule{C}) where {C <: RingElement}
    return [Symbol("v_$(i)") for i in 1:ngens(V)]
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

# no special ones


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_of_action(x::MatElem{C}, V::LieAlgebraStdModule{C}) where {C <: RingElement}
    size(x, 1) == size(x, 2) == V.n || throw(ArgumentError("Wrong size of matrix"))
    return x
end


###############################################################################
#
#   Constructor
#
###############################################################################

function standard_module(R::Ring, n::Int)
    return LieAlgebraStdModule{elem_type(R)}(R, n)
end
