function patch_oscar_serialization_namespace()
    Oscar.get_oscar_serialization_version() # call once to ensure the Oscar version is set

    push!(
        Oscar.oscar_serialization_version[],
        :PBWDeformations => ["https://github.com/lgoettgens/PBWDeformations.jl", Base.get_pkgversion_from_path(dirname(@__DIR__))],
    )
end

function register_serialization_types()
    # modifies Oscar global variables, so must be called from the __init__ function
    @eval @register_serialization_type SmashProductLie uses_id
    @eval @register_serialization_type SmashProductLieElem
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
