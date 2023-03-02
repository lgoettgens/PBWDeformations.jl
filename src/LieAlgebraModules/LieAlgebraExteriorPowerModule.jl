struct LieAlgebraExteriorPowerModule{C <: RingElement} <: LieAlgebraModule{C}
    inner_mod::LieAlgebraModule
    power::Int
    transformation_matrix_cache::Dict{MatElem{C}, MatElem{C}}
    ind_map::Vector{Vector{Int}}

    function LieAlgebraExteriorPowerModule{C}(inner_mod::LieAlgebraModule, power::Int) where {C <: RingElement}
        transformation_matrix_cache = Dict{MatElem{C}, MatElem{C}}()
        ind_map = collect(Combinatorics.combinations(1:ngens(inner_mod), power))
        return new{C}(inner_mod, power, transformation_matrix_cache, ind_map)
    end
end

struct LieAlgebraExteriorPowerModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraExteriorPowerModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraExteriorPowerModuleElem{C}}) where {C <: RingElement} = LieAlgebraExteriorPowerModule{C}

elem_type(::Type{LieAlgebraExteriorPowerModule{C}}) where {C <: RingElement} = LieAlgebraExteriorPowerModuleElem{C}

parent(v::LieAlgebraExteriorPowerModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement} = base_ring(V.inner_mod)

ngens(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement} = binomial(ngens(V.inner_mod), V.power)


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement}
    print(io, "$(V.power)-th tensor power of ")
    print(IOContext(io, :compact => true), V.inner_mod)
end

function symbols(V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement}
    if V.power == 1
        return symbols(V.inner_mod)
    end
    if isa(V.inner_mod, LieStdModule)
        parentheses = identity
    else
        parentheses = x -> "($x)"
    end

    return [join(s .|> parentheses, " ∧ ") for s in Combinatorics.combinations(symbols(V.inner_mod), V.power)]
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::LieAlgebraExteriorPowerModule{C})(
    a::Vector{T},
) where {T <: LieAlgebraModuleElem{C}} where {C <: RingElement}
    length(a) == V.power || error("Length of vector does not match tensor power.")
    all(x -> parent(x) == V.inner_mod, a) || error("Incompatible modules.")
    mat = zero_matrix(base_ring(V), 1, ngens(V))
    for (i, _inds) in enumerate(V.ind_map), inds in Combinatorics.permutations(_inds)
        sgn = levicivita(sortperm(inds))
        mat[1, i] += sgn * prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
    return LieAlgebraExteriorPowerModuleElem{C}(V, mat)
end


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_of_action(x::MatElem{C}, V::LieAlgebraExteriorPowerModule{C}) where {C <: RingElement}
    T = tensor_power(V.inner_mod, V.power)
    basis_change_E2T = zero(x, ngens(T), ngens(V))
    basis_change_T2E = zero(x, ngens(V), ngens(T))

    for (i, _inds) in enumerate(V.ind_map), inds in Combinatorics.permutations(_inds)
        sgn = levicivita(sortperm(inds))
        j = findfirst(==(inds), T.ind_map)
        basis_change_E2T[j, i] = sgn // factorial(V.power)
        basis_change_T2E[i, j] = sgn
    end

    xT = transformation_matrix_of_action(x, T)

    return basis_change_T2E * xT * basis_change_E2T
end


###############################################################################
#
#   Constructor
#
###############################################################################

function exterior_power(V::LieAlgebraModule{C}, k::Int) where {C <: RingElement}
    return LieAlgebraExteriorPowerModule{C}(V, k)
end
