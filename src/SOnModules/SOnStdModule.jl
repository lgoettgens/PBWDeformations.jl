struct SOnStdModule{C <: RingElement} <: SOnModule{C}
    R::Ring
    n::Int
end

struct SOnStdModuleElem{C <: RingElement} <: SOnModuleElem{C}
    parent::SOnStdModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{SOnStdModuleElem{C}}) where {C <: RingElement} = SOnStdModule{C}

elem_type(::Type{SOnStdModule{C}}) where {C <: RingElement} = SOnStdModuleElem{C}

parent(v::SOnStdModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::SOnStdModule{C}) where {C <: RingElement} = V.R

ngens(V::SOnStdModule{C}) where {C <: RingElement} = V.n


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::SOnStdModule{C}) where {C <: RingElement}
    print(io, "SO_$(V.n)_StdModule over ")
    print(IOContext(io, :compact => true), base_ring(V))
end

function symbols(V::SOnStdModule{C}) where {C <: RingElement}
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
#   Constructor
#
###############################################################################

function so_standard_module(R::Ring, n::Int)
    return SOnStdModule{elem_type(R)}(R, n)
end