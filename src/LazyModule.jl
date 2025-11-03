abstract type LazyLieAlgebraModule{C <: FieldElem, LieT <: LieAlgebraElem{C}} <: AbstractAlgebra.Set end

const LieAlgebraModuleOrLazy{C <: FieldElem, LieT <: LieAlgebraElem{C}} = Union{LieAlgebraModule{C, LieT}, LazyLieAlgebraModule{C, LieT}}

Oscar._is_standard_module(::LazyLieAlgebraModule) = false
Oscar._is_dual(::LazyLieAlgebraModule) = (false, nothing)
Oscar._is_direct_sum(::LazyLieAlgebraModule) = (false, nothing)
Oscar._is_tensor_product(::LazyLieAlgebraModule) = (false, nothing)
Oscar._is_exterior_power(::LazyLieAlgebraModule) = (false, nothing, nothing)
Oscar._is_symmetric_power(::LazyLieAlgebraModule) = (false, nothing, nothing)
Oscar._is_tensor_power(::LazyLieAlgebraModule) = (false, nothing, nothing)


################################################################################
#
# Tensor products
#
################################################################################

mutable struct LazyTensorProductLieAlgebraModule{
    C <: FieldElem,
    LieT <: LieAlgebraElem{C},
    FactorsT <: Tuple{Vararg{Union{LieAlgebraModule{C, LieT}, LazyLieAlgebraModule{C, LieT}}}},
} <: LazyLieAlgebraModule{C, LieT}
    Ls::FactorsT
end

function lazy_tensor_product(V::LieAlgebraModuleOrLazy{C, LieT}, Vs::LieAlgebraModuleOrLazy{C, LieT}...) where {C <: FieldElem, LieT <: LieAlgebraElem{C}}
    return LazyTensorProductLieAlgebraModule{C, LieT, typeof((V, Vs...))}((V, Vs...))
end

Oscar._is_tensor_product(L::LazyTensorProductLieAlgebraModule) = (true, L.Ls)


################################################################################
#
# Exterior powers
#
################################################################################

mutable struct LazyExteriorPowerLieAlgebraModule{
    C <: FieldElem,
    LieT <: LieAlgebraElem{C},
    InnerT <: Union{LieAlgebraModule{C, LieT}, LazyLieAlgebraModule{C, LieT}},
} <: LazyLieAlgebraModule{C, LieT}
    L::InnerT
    k::Int
end

function lazy_exterior_power_obj(L::LieAlgebraModuleOrLazy{C, LieT}, k::Int) where {C <: FieldElem, LieT <: LieAlgebraElem{C}}
    return LazyExteriorPowerLieAlgebraModule{C, LieT, typeof(L)}(L, k)
end

Oscar._is_exterior_power(L::LazyExteriorPowerLieAlgebraModule) = (true, L.L, L.k)
