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
    setup_filepath = joinpath(path, string_for_filename_setup(sp)) * file_ext
    if isfile(setup_filepath)
        @vprint :PBWDeformationsDatabase "Setup file already exists, loading it..."
        sp, _ = load(setup_filepath)::Tuple{typeof(sp), MatSpace{elem_type(sp)}}
        @vprintln :PBWDeformationsDatabase " Done"
    else
        @vprint :PBWDeformationsDatabase "Saving setup file..."
        save(joinpath(path, string_for_filename_setup(sp)) * file_ext, (sp, parent(zero_matrix(sp, dim(base_module(sp)), dim(base_module(sp))))))
        @vprintln :PBWDeformationsDatabase " Done"
    end
    for deg in 0:maxdeg
        for degs in [deg:deg, 0:deg]
            @vprintln :PBWDeformationsDatabase "Computing GlnGraphDeformBasis for degrees $(degs)..." # ProgressMeter needs a linebreak
            b = GlnGraphDeformBasis(sp, degs)
            @vprintln :PBWDeformationsDatabase " Done"
            @vprint :PBWDeformationsDatabase "Saving GlnGraphDeformBasis..."
            save(joinpath(path, string_for_filename(b)) * file_ext, b; serializer=PBWDeformations.JSONSerializerNoRefs())
            @vprintln :PBWDeformationsDatabase " Done"
            @vprint :PBWDeformationsDatabase "Computing PBW deformations for degrees $(degs)..."
            ms = all_pbwdeformations(sp, b)
            @vprintln :PBWDeformationsDatabase " Done"
            @vprint :PBWDeformationsDatabase "Saving PBW deformations..."
            save(joinpath(path, string_for_filename_pbwdeforms(b)) * file_ext, ms; serializer=PBWDeformations.JSONSerializerNoRefs())
            @vprintln :PBWDeformationsDatabase " Done"
        end
    end
end

function prepare_loading(base_path::String, sp::SmashProductLie)
    @req isdir(base_path) "Base path must be an existing directory"
    path = joinpath(base_path, string_for_path(sp))
    @req isdir(path) "This instance does not exist in the database"
    setup_filepath = joinpath(path, string_for_filename_setup(sp)) * file_ext
    @req isfile(setup_filepath) "This instance does not exist in the database"
    @vprint :PBWDeformationsDatabase "Setup file exists, loading it..."
    sp, _ = load(setup_filepath)::Tuple{typeof(sp), MatSpace{elem_type(sp)}}
    @vprintln :PBWDeformationsDatabase " Done"
    return path
end

function load_glngraph_deform_basis(base_path::String, sp::SmashProductLie, degs::AbstractVector{Int})
    path = prepare_loading(base_path, sp)

    filepath = joinpath(path, string_for_filename(GlnGraphDeformBasis, sp, degs)) * file_ext
    @req isfile(filepath) "The requested degree does not exist in the database"
    @vprint :PBWDeformationsDatabase "Found GlnGraphDeformBasis for degree $(degs). Loading..."
    b = load(filepath)::GlnGraphDeformBasis{elem_type(sp)}
    @vprintln :PBWDeformationsDatabase " Done"
    return b
end

function load_glngraph_deform_bases(base_path::String, sp::SmashProductLie, degss::AbstractVector{<:AbstractVector{Int}})
    path = prepare_loading(base_path, sp)

    bs = GlnGraphDeformBasis{elem_type(sp)}[]
    for degs in degss
        filepath = joinpath(path, string_for_filename(GlnGraphDeformBasis, sp, degs)) * file_ext
        @req isfile(filepath) "The requested degree does not exist in the database"
        @vprint :PBWDeformationsDatabase "Found GlnGraphDeformBasis for degree $(degs). Loading..."
        push!(bs, load(filepath)::GlnGraphDeformBasis{elem_type(sp)})
        @vprintln :PBWDeformationsDatabase " Done"
    end
    return bs
end

function load_glngraph_deform_bases(base_path::String, sp::SmashProductLie)
    path = prepare_loading(base_path, sp)

    deg = 0
    bs = GlnGraphDeformBasis{elem_type(sp)}[]
    while (degs = deg:deg; filepath = joinpath(path, string_for_filename(GlnGraphDeformBasis, sp, degs)) * file_ext; isfile(filepath))
        @vprint :PBWDeformationsDatabase "Found GlnGraphDeformBasis for degree $(degs). Loading..."
        push!(bs, load(filepath)::GlnGraphDeformBasis{elem_type(sp)})
        @vprintln :PBWDeformationsDatabase " Done"
        deg += 1
    end
    return bs
end

function load_pbwdeformations(base_path::String, sp::SmashProductLie, degs::AbstractVector{Int})
    path = prepare_loading(base_path, sp)

    filepath = joinpath(path, string_for_filename_pbwdeforms(sp, degs)) * file_ext
    @req isfile(filepath) "The requested degree does not exist in the database"
    @vprint :PBWDeformationsDatabase "Found PBW deformations for degree $(degs). Loading..."
    ms = load(filepath)::Vector{DeformationMap{elem_type(sp)}}
    @vprintln :PBWDeformationsDatabase " Done"
    return ms
end

function load_pbwdeformations(base_path::String, sp::SmashProductLie, degss::AbstractVector{<:AbstractVector{Int}})
    path = prepare_loading(base_path, sp)

    ms = Vector{DeformationMap{elem_type(sp)}}[]
    for degs in degss
        filepath = joinpath(path, string_for_filename_pbwdeforms(sp, degs)) * file_ext
        @req isfile(filepath) "The requested degree does not exist in the database"
        @vprint :PBWDeformationsDatabase "Found PBW deformations for degree $(degs). Loading..."
        push!(ms, load(filepath)::Vector{DeformationMap{elem_type(sp)}})
        @vprintln :PBWDeformationsDatabase " Done"
    end
    return ms
end

function load_pbwdeformations(base_path::String, sp::SmashProductLie)
    path = prepare_loading(base_path, sp)

    deg = 0
    mss = Vector{DeformationMap{elem_type(sp)}}[]
    while (degs = deg:deg; filepath = joinpath(path, string_for_filename_pbwdeforms(sp, degs)) * file_ext; isfile(filepath))
        @vprint :PBWDeformationsDatabase "Found PBW deformations for degree $(degs). Loading..."
        push!(mss, Vector{DeformationMap{elem_type(sp)}}(load(filepath))::Vector{DeformationMap{elem_type(sp)}}) # see https://github.com/oscar-system/Oscar.jl/issues/3983
        @vprintln :PBWDeformationsDatabase " Done"
        deg += 1
    end
    return mss
end


################################################################################
#
# File path helpers
#
################################################################################

const file_ext = ".mrdi"

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
            return "DV"
        end
        return "D_" * string_for_filename(B) * "_"
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
    return string_for_filename(typeof(b), b.sp, b.degs)
end

function string_for_filename(::Type{<:ArcDiagDeformBasis}, sp::SmashProductLie, degs::AbstractVector{Int})
    return "ArcDiagDeformBasis-" * string_for_filename(sp) * "-" * join(degs, "_")
end

function string_for_filename(b::GlnGraphDeformBasis)
    return string_for_filename(typeof(b), b.sp, b.degs)
end

function string_for_filename(::Type{<:GlnGraphDeformBasis}, sp::SmashProductLie, degs::AbstractVector{Int})
    return "GlnGraphDeformBasis-" * string_for_filename(sp) * "-" * join(degs, "_")
end


function string_for_filename_pbwdeforms(b::Union{ArcDiagDeformBasis, GlnGraphDeformBasis})
    return string_for_filename_pbwdeforms(b.sp, b.degs)
end

function string_for_filename_pbwdeforms(sp::SmashProductLie, degs::AbstractVector{Int})
    return "PBWDeformations-" * string_for_filename(sp) * "-" * join(degs, "_")
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
