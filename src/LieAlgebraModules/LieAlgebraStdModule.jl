struct LieAlgebraStdModule{C <: RingElement} <: LieAlgebraModule{C}
    L::LieAlgebra{C}

    function LieAlgebraStdModule{C}(L::LieAlgebra{C}) where {C <: RingElement}
        return new{C}(L)
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

base_ring(V::LieAlgebraStdModule{C}) where {C <: RingElement} = base_ring(base_liealgebra(V))

base_liealgebra(V::LieAlgebraStdModule{C}) where {C <: RingElement} = V.L

ngens(V::LieAlgebraStdModule{C}) where {C <: RingElement} = base_liealgebra(V).n


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraStdModule{C}) where {C <: RingElement}
    print(io, "StdModule of ")
    print(IOContext(io, :compact => true), base_liealgebra(V))
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
#   Comparison functions
#
###############################################################################

function Base.:(==)(V1::LieAlgebraStdModule{C}, V2::LieAlgebraStdModule{C}) where {C <: RingElement}
    return (V1.L) == (V2.L)
end


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_by_basisindex(V::LieAlgebraStdModule{C}, i::Int) where {C <: RingElement}
    return basis(base_liealgebra(V))[i]
end


###############################################################################
#
#   Constructor
#
###############################################################################

function standard_module(L::LieAlgebra{C}) where {C <: RingElement}
    return LieAlgebraStdModule{elem_type(base_ring(L))}(L)
end
