struct SOnTensorPowerModule{C <: RingElement} <: SOnModule{C}
    inner_mod::SOnModule
    power::Int
    ind_map::Vector{Vector{Int}}

    function SOnTensorPowerModule{C}(inner_mod::SOnModule, power::Int) where {C <: RingElement}
        ind_map = collect(ProductIterator(1:ngens(inner_mod), power))
        return new{C}(inner_mod, power, ind_map)
    end
end

struct SOnTensorPowerModuleElem{C <: RingElement} <: SOnModuleElem{C}
    parent::SOnTensorPowerModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{SOnTensorPowerModuleElem{C}}) where {C <: RingElement} = SOnTensorPowerModule{C}

elem_type(::Type{SOnTensorPowerModule{C}}) where {C <: RingElement} = SOnTensorPowerModuleElem{C}

parent(v::SOnTensorPowerModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::SOnTensorPowerModule{C}) where {C <: RingElement} = base_ring(V.inner_mod)

ngens(V::SOnTensorPowerModule{C}) where {C <: RingElement} = ngens(V.inner_mod)^V.power


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::SOnTensorPowerModule{C}) where {C <: RingElement}
    print(io, "$(V.power)-th tensor power of ")
    print(IOContext(io, :compact => true), V.inner_mod)
end

function symbols(V::SOnTensorPowerModule{C}) where {C <: RingElement}
    if V.power == 1
        return symbols(V.inner_mod)
    end
    if isa(V.inner_mod, SOnStdModule)
        parentheses = identity
    else
        parentheses = x -> "($x)"
    end

    return [join(s .|> parentheses, " ⊗ ") for s in ProductIterator(symbols(V.inner_mod), V.power)]
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::SOnTensorPowerModule{C})(a::Vector{T}) where {T <: SOnModuleElem{C}} where {C <: RingElement}
    length(a) == V.power || error("Length of vector does not match tensor power.")
    all(x -> parent(x) == V.inner_mod, a) || error("Incompatible modules.")
    mat = zero_matrix(base_ring(V), 1, ngens(V))
    for (i, es) in enumerate(ProductIterator([a[i].mat for i in 1:length(a)]))
        mat[1, i] = prod(es)
    end
    return SOnTensorPowerModuleElem{C}(V, mat)
end


###############################################################################
#
#   Constructor
#
###############################################################################

function tensor_power(V::SOnModule{C}, k::Int) where {C <: RingElement}
    return SOnTensorPowerModule{C}(V, k)
end
