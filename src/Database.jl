module Database

using ..PBWDeformations
using Oscar

using Oscar: _is_dual
using Oscar: _is_direct_sum
using Oscar: _is_exterior_power
using Oscar: _is_symmetric_power
using Oscar: _is_tensor_power
using Oscar: _is_tensor_product
using Oscar: _is_standard_module

################################################################################
#
# Main function
#
################################################################################

function compute_and_save_instance(base_path::String, sp::SmashProductLie, maxdeg::Int)
    path = joinpath(base_path, string_for_path(sp))
    mkpath(path)
    setup_filepath = joinpath(path, string_for_filename_setup(sp)) * ".mrdi"
    if isfile(setup_filepath)
        @vprintln :PBWDeformationsDatabase "Setup file already exists, loading it..."
        sp, _ = load(setup_filepath)::Tuple{typeof(sp), MatSpace{elem_type(sp)}}
        @vprintln :PBWDeformationsDatabase " Done"
    else
        @vprintln :PBWDeformationsDatabase "Saving setup file..."
        save(joinpath(path, string_for_filename_setup(sp))* ".mrdi", (sp, parent(zero_matrix(sp, dim(base_module(sp)), dim(base_module(sp))))))
        @vprintln :PBWDeformationsDatabase " Done"
    end
    for deg in 0:maxdeg
        for degs in [deg:deg, 0:deg]
            @vprintln :PBWDeformationsDatabase "Computing GlnGraphDeformBasis for degrees $(degs)..."
            b = GlnGraphDeformBasis(sp, degs)
            @vprintln :PBWDeformationsDatabase " Done"
            @vprintln :PBWDeformationsDatabase "Saving GlnGraphDeformBasis..."
            save(joinpath(path, string_for_filename(b)) * ".mrdi", b; serializer=PBWDeformations.JSONSerializerNoRefs())
            @vprintln :PBWDeformationsDatabase " Done"
            @vprintln :PBWDeformationsDatabase "Computing PBW deformations for degrees $(degs)..."
            ms = all_pbwdeformations(sp, b)
            @vprintln :PBWDeformationsDatabase " Done"
            @vprintln :PBWDeformationsDatabase "Saving PBW deformations..."
            save(joinpath(path, string_for_filename_pbwdeforms(b))* ".mrdi", ms; serializer=PBWDeformations.JSONSerializerNoRefs())
            @vprintln :PBWDeformationsDatabase " Done"
        end
    end
end


################################################################################
#
# File path helpers
#
################################################################################

function type_as_string(L::LinearLieAlgebra)
    type = get_attribute(L, :type, nothing)
    if type == :general_linear
        return "gl"
    elseif type == :special_orthogonal
        return "so"
    end
    error("type_as_string not implemented for this type of Lie algebra")
end

function string_for_filename(L::LinearLieAlgebra)
    type = type_as_string(L)
    n = L.n
    return "$(type)_$(n)"
end

function string_for_filename(V::LieAlgebraModule)
    if _is_standard_module(V)
        return "V"
    elseif ((fl, B) = _is_dual(V); fl)
        if _is_standard_module(B)
            return "Vd"
        end
        return "dual_" * string_for_filename(B) * "_"
    elseif ((fl, Bs) = _is_direct_sum(V); fl)
        return "_" * join(string_for_filename.(Bs), "_+_") * "_"
    elseif ((fl, Bs) = _is_tensor_product(V); fl)
        return "_" * join(string_for_filename.(Bs), "_x_") * "_"
    elseif ((fl, B, k) = _is_exterior_power(V); fl)
        return "E$(k)_" * string_for_filename(B) * "_"
    elseif ((fl, B, k) = _is_symmetric_power(V); fl)
        return "S$(k)_" * string_for_filename(B) * "_"
    elseif ((fl, B, k) = _is_tensor_power(V); fl)
        return "T$(k)_" * string_for_filename(B) * "_"
    end
    error("not implemented for this type of module")
end

function string_for_filename(sp::SmashProductLie)
    return string_for_filename(base_lie_algebra(sp)) * "-" * string_for_filename(base_module(sp))
end

function string_for_filename(b::ArcDiagDeformBasis)
    return "ArcDiagDeformBasis-" * string_for_filename(b.sp) * join(b.degs, "_")
end

function string_for_filename(b::GlnGraphDeformBasis)
    return "GlnGraphDeformBasis-" * string_for_filename(b.sp) * join(b.degs, "_")
end


function string_for_filename_pbwdeforms(b::Union{ArcDiagDeformBasis, GlnGraphDeformBasis})
    return "PBWDeformations-" * string_for_filename(b.sp) * join(b.degs, "_")
end


function string_for_filename_setup(b::Union{ArcDiagDeformBasis, GlnGraphDeformBasis})
    return string_for_filename_setup(b.sp)
end

function string_for_filename_setup(sp::SmashProductLie)
    return "setup-" * string_for_filename(sp)
end


function string_for_path(b::Union{ArcDiagDeformBasis, GlnGraphDeformBasis})
    return string_for_path(b.sp)
end

function string_for_path(sp::SmashProductLie)
    L = base_lie_algebra(sp)
    V = base_module(sp)
    return joinpath(type_as_string(L), string_for_filename(V))
end

end # module
