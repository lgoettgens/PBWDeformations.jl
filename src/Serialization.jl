function patch_oscar_serialization_namespace()
    Oscar.get_oscar_serialization_version() # call once to ensure the Oscar version is set

    push!(
        Oscar.oscar_serialization_version[],
        :PBWDeformations => ["https://github.com/lgoettgens/PBWDeformations.jl", Base.get_pkgversion_from_path(dirname(@__DIR__))],
    )
end

function register_serialization_types()
    @eval @register_serialization_type GlnGraph

    @eval @register_serialization_type SmashProductLie uses_id
    @eval @register_serialization_type SmashProductLieElem

    @eval @register_serialization_type SmashProductLieDeform uses_id
    @eval @register_serialization_type SmashProductLieDeformElem
end

###############################################################################
#
#   GlnGraph
#
###############################################################################

type_params(g::GlnGraph) = TypeParams(GlnGraph, nothing)


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
        save_object(s, e.simplified, :is_simplified)
    end
end

function load_object(s::DeserializerState, ::Type{<:SmashProductLieElem}, sp::SmashProductLie)
    e = sp(load_object(s, elem_type(underlying_algebra(sp)), underlying_algebra(sp), :data))
    e.simplified = load_object(s, Bool, :is_simplified)
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
        save_object(s, e.simplified, :is_simplified)
    end
end

function load_object(s::DeserializerState, ::Type{<:SmashProductLieDeformElem}, d::SmashProductLieDeform)
    e = d(load_object(s, elem_type(underlying_algebra(d)), underlying_algebra(d), :data))
    e.simplified = load_object(s, Bool, :is_simplified)
    return e
end
