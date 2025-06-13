const ArcDiagDeformBasisDataT = ArcDiagram

"""
Concrete subtype of [`DeformBasis`](@ref).
Each element of the basis is induced by an arc diagram of a suitable size,
which gets symmetrized and specialised to the given smash product.
This process is due to [FM22](@cite).
"""
const ArcDiagDeformBasis{T} = ArcDiagBasedDeformBasis{ArcDiagDeformBasisDataT, T} where {T <: SmashProductLieElem}

function ArcDiagDeformBasis(
    sp::SmashProductLie{C, LieC, LieT},
    degs::AbstractVector{Int};
    no_normalize::Bool=false,
) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    LieType = Val(get_attribute(base_lie_algebra(sp), :type, nothing)::Union{Nothing, Symbol})
    @req LieType isa Union{SO, GL} "Only works for so_n and gl_n."
    if LieType isa SO && has_attribute(base_lie_algebra(sp), :form)
        @req isone(get_attribute(base_lie_algebra(sp), :form)::dense_matrix_type(C)) "Only works for so_n represented as skew-symmetric matrices."
    end
    return ArcDiagDeformBasis(LieType, sp, degs; no_normalize)
end

function data_iter_and_len(
    ::Type{ArcDiagDeformBasis},
    LieType::Union{SO, GL},
    W::LieAlgebraModule,
    case::Symbol,
    d::Int,
)
    diag_iter = pbw_arc_diagrams(LieType, W, d)
    len = length(diag_iter)
    return diag_iter, len::Int
end

function should_data_be_used(
    ::Type{ArcDiagDeformBasis},
    LieType::Union{SO, GL},
    data::ArcDiagDeformBasisDataT,
    ::SmashProductLie,
    ::LieAlgebraModule,
    ::Symbol,
    cache::Union{Dict{<:Any, Bool}, Nothing},
)
    diag = data
    is_crossing_free(diag; part=:lower)
end


function pbw_arc_diagrams(T::Union{SO, GL}, V::LieAlgebraModule, d::Int)
    upper_verts = arc_diagram_upper_points(T, V)
    lower_verts = arc_diagram_lower_points(T, V, d)
    upper_iss = arc_diagram_upper_iss(T, V)
    lower_iss = arc_diagram_lower_iss(T, V, d)
    indep_sets = Vector{Int}[[(-1) .* is for is in upper_iss]; [is for is in lower_iss]]
    return all_arc_diagrams(arc_diagram_type(T), upper_verts, lower_verts; indep_sets)
end


arc_diagram_type(::SO) = Undirected

arc_diagram_type(::GL) = Directed


function is_tensor_generator(V::LieAlgebraModule)
    if _is_standard_module(V)
        return true
    end
    fl, base = _is_dual(V)
    return fl && _is_standard_module(base)
end


function arc_diagram_upper_points(T::SO, V::LieAlgebraModule)
    if _is_standard_module(V)
        return 1
    elseif ((fl, Ws) = _is_tensor_product(V); fl)
        return sum(arc_diagram_upper_points(T, W) for W in Ws)
    elseif ((fl, W, k) = is_power_with_data(V); fl)
        return arc_diagram_upper_points(T, W) * k
    else
        error("Not implemented.")
    end
end

function arc_diagram_upper_points(T::GL, V::LieAlgebraModule)
    if _is_standard_module(V)
        return [true]
    elseif ((fl, W) = _is_dual(V); fl) && _is_standard_module(W)
        return [false]
    elseif ((fl, Ws) = _is_tensor_product(V); fl)
        return reduce(vcat, [arc_diagram_upper_points(T, W) for W in Ws])
    elseif ((fl, W, k) = is_power_with_data(V); fl)
        upper_points_W = arc_diagram_upper_points(T, W)
        return reduce(vcat, [upper_points_W for _ in 1:k])
    else
        error("Not implemented.")
    end
end

function arc_diagram_num_upper_points(T::SO, V::LieAlgebraModule)
    return arc_diagram_upper_points(T, V)
end

function arc_diagram_num_upper_points(T::GL, V::LieAlgebraModule)
    return length(arc_diagram_upper_points(T, V))
end


function arc_diagram_upper_iss(T::Union{SO, GL}, V::LieAlgebraModule)
    if is_tensor_generator(V)
        return Vector{Int}[]
    elseif ((fl, inner_mods) = _is_tensor_product(V); fl)
        offset = 0
        iss = Vector{Int}[]
        for mod in inner_mods
            append!(iss, [is .+ offset for is in arc_diagram_upper_iss(T, mod)])
            offset += arc_diagram_num_upper_points(T, mod)
        end
        return iss
    elseif ((fl, inner_mod, power) = is_power_with_data(V); fl)
        if is_tensor_generator(inner_mod)
            if _is_exterior_power(V)[1]
                return [collect(1:power)]
            else
                return Vector{Int}[]
            end
        else
            iss = arc_diagram_upper_iss(T, inner_mod)
            return [is .+ k * arc_diagram_num_upper_points(T, inner_mod) for k in 0:power-1 for is in iss]
        end
    else
        error("Not implemented.")
    end
end


function arc_diagram_lower_points(::SO, _::LieAlgebraModule, d::Int)
    # L ≅ Sᵈ ⋀² V
    return 2d
end

function arc_diagram_lower_points(::GL, _::LieAlgebraModule, d::Int)
    # L ≅ Sᵈ (V ⊗ V*)
    return reduce(vcat, ([1, 0] for _ in 1:d); init=Int[])
end

function arc_diagram_num_lower_points(T::SO, V::LieAlgebraModule, d::Int)
    return arc_diagram_lower_points(T, V, d)
end

function arc_diagram_num_lower_points(T::GL, V::LieAlgebraModule, d::Int)
    return length(arc_diagram_lower_points(T, V, d))
end


function arc_diagram_lower_iss(::SO, _::LieAlgebraModule, d::Int)
    # L ≅ Sᵈ ⋀² V
    return collect([2i - 1, 2i] for i in 1:d)
end

function arc_diagram_lower_iss(::GL, _::LieAlgebraModule, _::Int)
    # L ≅ Sᵈ (V ⊗ V*)
    return Vector{Int}[]
end
