struct LieAlgebraSymmetricPowerModule{C <: RingElement} <: LieAlgebraModule{C}
    inner_mod::LieAlgebraModule
    power::Int
    transformation_matrix_cache::Dict{LieAlgebraElem{C}, MatElem{C}}
    ind_map::Vector{Vector{Int}}

    function LieAlgebraSymmetricPowerModule{C}(inner_mod::LieAlgebraModule{C}, power::Int) where {C <: RingElement}
        transformation_matrix_cache = Dict{LieAlgebraElem{C}, MatElem{C}}()
        ind_map = collect(Combinatorics.with_replacement_combinations(1:ngens(inner_mod), power))
        return new{C}(inner_mod, power, transformation_matrix_cache, ind_map)
    end
end

struct LieAlgebraSymmetricPowerModuleElem{C <: RingElement} <: LieAlgebraModuleElem{C}
    parent::LieAlgebraSymmetricPowerModule{C}
    mat::MatElem{C}
end


###############################################################################
#
#   Basic manipulation
#
###############################################################################

parent_type(::Type{LieAlgebraSymmetricPowerModuleElem{C}}) where {C <: RingElement} = LieAlgebraSymmetricPowerModule{C}

elem_type(::Type{LieAlgebraSymmetricPowerModule{C}}) where {C <: RingElement} = LieAlgebraSymmetricPowerModuleElem{C}

parent(v::LieAlgebraSymmetricPowerModuleElem{C}) where {C <: RingElement} = v.parent

base_ring(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement} = base_ring(V.inner_mod)

ngens(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement} =
    binomial(ngens(V.inner_mod) + V.power - 1, V.power)


###############################################################################
#
#   String I/O
#
###############################################################################

function Base.show(io::IO, V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement}
    print(io, "$(V.power)-th tensor power of ")
    print(IOContext(io, :compact => true), V.inner_mod)
end

function symbols(V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement}
    if V.power == 1
        return symbols(V.inner_mod)
    end
    if isa(V.inner_mod, LieStdModule)
        parentheses = identity
    else
        parentheses = x -> "($x)"
    end

    return [
        begin
            join((
                begin
                    e = count(==(i), inds)
                    if e == 1
                        s |> parentheses
                    else
                        "$(s |> parentheses)^$e"
                    end
                end for (i, s) in enumerate(symbols(V.inner_mod)) if in(i, inds)
            ), "*")
        end for inds in Combinatorics.with_replacement_combinations(1:ngens(V.inner_mod), V.power)
    ]
end


###############################################################################
#
#   Parent object call overload
#
###############################################################################

function (V::LieAlgebraSymmetricPowerModule{C})(
    a::Vector{T},
) where {T <: LieAlgebraModuleElem{C}} where {C <: RingElement}
    length(a) == V.power || error("Length of vector does not match tensor power.")
    all(x -> parent(x) == V.inner_mod, a) || error("Incompatible modules.")
    mat = zero_matrix(base_ring(V), 1, ngens(V))
    for (i, _inds) in enumerate(V.ind_map), inds in unique(Combinatorics.permutations(_inds))
        mat[1, i] += prod(a[j].mat[k] for (j, k) in enumerate(inds))
    end
    return LieAlgebraSymmetricPowerModuleElem{C}(V, mat)
end


###############################################################################
#
#   Comparison functions
#
###############################################################################

function Base.:(==)(
    V1::LieAlgebraSymmetricPowerModule{C},
    V2::LieAlgebraSymmetricPowerModule{C},
) where {C <: RingElement}
    return (V1.inner_mod, V1.power) == (V2.inner_mod, V2.power)
end


###############################################################################
#
#   Module action
#
###############################################################################

function transformation_matrix_of_action(x::MatElem{C}, V::LieAlgebraSymmetricPowerModule{C}) where {C <: RingElement}
    T = tensor_power(V.inner_mod, V.power)
    basis_change_S2T = zero(x, ngens(T), ngens(V))
    basis_change_T2S = zero(x, ngens(V), ngens(T))

    for (i, _inds) in enumerate(V.ind_map), inds in Combinatorics.permutations(_inds)
        j = findfirst(==(inds), T.ind_map)
        basis_change_S2T[j, i] += 1 // factorial(V.power)
        basis_change_T2S[i, j] = 1
    end

    xT = transformation_matrix_of_action(x, T)

    return basis_change_T2S * xT * basis_change_S2T
end


###############################################################################
#
#   Constructor
#
###############################################################################

function symmetric_power(V::LieAlgebraModule{C}, k::Int) where {C <: RingElement}
    return LieAlgebraSymmetricPowerModule{C}(V, k)
end
