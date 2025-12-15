@static if !isdefined(Oscar, :Serialization) # introduced in https://github.com/oscar-system/Oscar.jl/pull/5024 in Oscar v1.5
    function patch_oscar_serialization_namespace()
        Oscar.get_oscar_serialization_version() # call once to ensure the Oscar version is set

        push!(
            Oscar.oscar_serialization_version[],
            :PBWDeformations => ["https://github.com/lgoettgens/PBWDeformations.jl", VERSION_NUMBER],
        )
    end
else
    function patch_oscar_serialization_namespace()
        Oscar.Serialization.get_oscar_serialization_version() # call once to ensure the Oscar version is set

        push!(
            Oscar.Serialization.oscar_serialization_version[],
            :PBWDeformations => ["https://github.com/lgoettgens/PBWDeformations.jl", VERSION_NUMBER],
        )
    end
end

function register_serialization_types()
    @eval @register_serialization_type Partition

    @eval @register_serialization_type ArcDiagramUndirected
    @eval @register_serialization_type ArcDiagramDirected

    @eval @register_serialization_type GlnGraph

    @eval @register_serialization_type SmashProductLie uses_id
    @eval @register_serialization_type SmashProductLieElem

    @eval @register_serialization_type SmashProductLieDeform uses_id
    @eval @register_serialization_type SmashProductLieDeformElem

    @eval @register_serialization_type StdDeformBasis

    @eval @register_serialization_type ArcDiagDeformBasis
    @eval @register_serialization_type GlnGraphDeformBasis
end

###############################################################################
#
#   ArcDiagram
#
###############################################################################

type_params(_::Partition{T}) where {T} = TypeParams(Partition, TypeParams(T, nothing))

function save_object(s::SerializerState, p::Partition)
    save_object(s, data(p))
end

function load_object(s::DeserializerState, ::Type{<:Partition{T}}) where {T}
    return partition(load_object(s, Vector{T}))
end

###############################################################################
#
#   ArcDiagram
#
###############################################################################

type_params(_::ArcDiagramUndirected) = TypeParams(ArcDiagramUndirected, nothing)

function save_object(s::SerializerState, a::ArcDiagramUndirected)
    save_data_dict(s) do
        save_object(s, n_upper_vertices(a), :n_upper_verts)
        save_object(s, n_lower_vertices(a), :n_lower_verts)
        save_object(s, a.upper_neighbors, :upper_neighbors)
        save_object(s, a.lower_neighbors, :lower_neighbors)
    end
end

function load_object(s::DeserializerState, ::Type{<:ArcDiagramUndirected})
    n_upper_verts = load_object(s, Int, :n_upper_verts)
    n_lower_verts = load_object(s, Int, :n_lower_verts)
    upper_neighbors = load_object(s, Vector{ArcDiagramVertex}, :upper_neighbors)
    lower_neighbors = load_object(s, Vector{ArcDiagramVertex}, :lower_neighbors)
    return arc_diagram(Undirected, n_upper_verts, n_lower_verts, upper_neighbors, lower_neighbors; check=false)
end

type_params(_::ArcDiagramDirected) = TypeParams(ArcDiagramDirected, nothing)

function save_object(s::SerializerState, a::ArcDiagramDirected)
    save_data_dict(s) do
        save_object(s, n_upper_vertices(a), :n_upper_verts)
        save_object(s, n_lower_vertices(a), :n_lower_verts)
        save_object(s, a.parity_upper_verts, :parity_upper_verts)
        save_object(s, a.parity_lower_verts, :parity_lower_verts)
        save_object(s, a.upper_neighbors, :upper_neighbors)
        save_object(s, a.lower_neighbors, :lower_neighbors)
    end
end

function load_object(s::DeserializerState, ::Type{<:ArcDiagramDirected})
    n_upper_verts = load_object(s, Int, :n_upper_verts)
    n_lower_verts = load_object(s, Int, :n_lower_verts)
    parity_upper_verts = load_object(s, Vector{Bool}, :parity_upper_verts)
    parity_lower_verts = load_object(s, Vector{Bool}, :parity_lower_verts)
    upper_neighbors = load_object(s, Vector{ArcDiagramVertex}, :upper_neighbors)
    lower_neighbors = load_object(s, Vector{ArcDiagramVertex}, :lower_neighbors)
    return arc_diagram(Directed, n_upper_verts, n_lower_verts, parity_upper_verts, parity_lower_verts,
                       upper_neighbors, lower_neighbors; check=false)
end

function load_object(s::DeserializerState, ::Type{ArcDiagram})
    if haskey(s, :parity_upper_verts)
        return load_object(s, ArcDiagramDirected)
    else
        return load_object(s, ArcDiagramUndirected)
    end
end

###############################################################################
#
#   GlnGraph
#
###############################################################################

type_params(_::GlnGraph) = TypeParams(GlnGraph, nothing)

function save_object(s::SerializerState, g::GlnGraph)
    save_data_dict(s) do
        save_object(s, g.n_left_verts, :n_left_verts)
        save_object(s, g.n_right_verts, :n_right_verts)
        save_object(s, g.parity_verts, :parity_verts)
        save_object(s, g.edges, :edges)
    end
end

function load_object(s::DeserializerState, ::Type{<:GlnGraph})
    n_left_verts = load_object(s, Int, :n_left_verts)
    n_right_verts = load_object(s, Int, :n_right_verts)
    parity_verts = load_object(s, Vector{Bool}, :parity_verts)
    edges = load_object(s, Vector{GlnGraphEdge}, :edges)
    return GlnGraph(n_left_verts, n_right_verts, parity_verts, edges; check=false, sort=false)
end


###############################################################################
#
#   SmashProductLie
#
###############################################################################

function type_params(sp::SmashProductLie)
    return TypeParams(
        SmashProductLie,
        :coefficient_ring => coefficient_ring(sp),
        :base_lie_algebra => base_lie_algebra(sp),
        :base_module => base_module(sp),
    )
end

function save_object(s::SerializerState, sp::SmashProductLie)
    save_data_dict(s) do
        # no data required but we leave this function here to generate a valid json
        # and leave root for possible future attrs
    end
end

function load_object(s::DeserializerState, ::Type{<:SmashProductLie}, d::Dict)
    R = d[:coefficient_ring]
    L = d[:base_lie_algebra]
    V = d[:base_module]
    return smash_product(R, L, V)
end

type_params(e::SmashProductLieElem) = TypeParams(SmashProductLieElem, parent(e))

function save_object(s::SerializerState, e::SmashProductLieElem)
    save_data_dict(s) do
        save_object(s, data(e), :data)
        if !is_zero(data(e))
            save_object(s, e.simplified, :is_simplified)
        end
    end
end

function load_object(s::DeserializerState, ::Type{<:SmashProductLieElem}, sp::SmashProductLie)
    e = sp(load_object(s, elem_type(underlying_algebra(sp)), underlying_algebra(sp), :data))
    if haskey(s, :is_simplified)
        e.simplified = load_object(s, Bool, :is_simplified)
    end
    return e
end


###############################################################################
#
#   SmashProductLieDeform
#
###############################################################################

function type_params(d::SmashProductLieDeform)
    return TypeParams(SmashProductLieDeform, d.sp)
end

function save_object(s::SerializerState, d::SmashProductLieDeform)
    save_data_dict(s) do
        save_object(s, d.kappa, :kappa)
    end
end

function load_object(s::DeserializerState, ::Type{<:SmashProductLieDeform}, sp::SmashProductLie{C, LieC, LieT}) where {C <: RingElem, LieC <: FieldElem, LieT <: LieAlgebraElem{LieC}}
    mat_space = matrix_space(sp, dim(base_module(sp)), dim(base_module(sp)))
    kappa = load_object(s, DeformationMap{SmashProductLieElem{C, LieC, LieT}}, mat_space, :kappa)
    return deform(sp, kappa)
end

type_params(e::SmashProductLieDeformElem) = TypeParams(SmashProductLieDeformElem, parent(e))

function save_object(s::SerializerState, e::SmashProductLieDeformElem)
    save_data_dict(s) do
        save_object(s, data(e), :data)
        if !is_zero(data(e))
            save_object(s, e.simplified, :is_simplified)
        end
    end
end

function load_object(s::DeserializerState, ::Type{<:SmashProductLieDeformElem}, d::SmashProductLieDeform)
    e = d(load_object(s, elem_type(underlying_algebra(d)), underlying_algebra(d), :data))
    if haskey(s, :is_simplified)
        e.simplified = load_object(s, Bool, :is_simplified)
    end
    return e
end

###############################################################################
#
#   StdDeformBasis
#
###############################################################################

function type_params(b::StdDeformBasis)
    return TypeParams(StdDeformBasis, b.sp)
end

function save_object(s::SerializerState, b::StdDeformBasis)
    save_data_dict(s) do
        save_object(s, b.degs, :degs)
    end
end

function load_object(s::DeserializerState, ::Type{<:StdDeformBasis}, sp::SmashProductLie)
    degs = load_object(s, Vector{Int}, :degs)
    return StdDeformBasis(sp, degs)
end

###############################################################################
#
#   ArcDiagBasedDeformBasis
#
###############################################################################

function type_params(b::ArcDiagBasedDeformBasis{ParamT}) where {ParamT}
    return TypeParams(ArcDiagBasedDeformBasis{ParamT}, b.sp)
end

function save_object(s::SerializerState, b::ArcDiagBasedDeformBasis)
    @req b.strict "Serialization is only supported for ArcDiagBasedDeformBasis that are constructor non-lazy"
    save_data_dict(s) do
        save_object(s, b.degs, :degs)
        save_object(s, length(b), :len)
        save_object(s, b.iter, :iter)
        save_object(s, b.param_reverse_map, :param_reverse_map)
        save_object(s, b.strict, :strict)
    end
end

function load_object(s::DeserializerState, ::Type{<:ArcDiagBasedDeformBasis{ParamT}}, sp::SmashProductLie) where {ParamT}
    degs = load_object(s, Vector{Int}, :degs)
    len = load_object(s, Int, :len)
    mat_space = matrix_space(sp, dim(base_module(sp)), dim(base_module(sp)))
    iter = load_object(s, Vector{DeformationMap{elem_type(sp)}}, mat_space, :iter)
    param_reverse_map = load_object(s, Dict{DeformationMap{elem_type(sp)}, Set{Tuple{Tuple{Int, Int}, ParamT}}}, Dict(:key_params => mat_space, :value_params => nothing), :param_reverse_map)
    strict = load_object(s, Bool, :strict)
    return ArcDiagBasedDeformBasis{ParamT, elem_type(sp)}(sp, degs, len, iter, param_reverse_map; strict)
end
