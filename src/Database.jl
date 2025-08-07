module Database

using ..PBWDeformations
using Oscar

using ..PBWDeformations: is_prefix_equal

using Oscar: _is_dual
using Oscar: _is_direct_sum
using Oscar: _is_exterior_power
using Oscar: _is_symmetric_power
using Oscar: _is_tensor_power
using Oscar: _is_tensor_product
using Oscar: _is_standard_module

################################################################################
#
# Auxiliary types
#
################################################################################

struct DBLieAlgebraKey
    type::Symbol
    n::Int
    F::String

    function DBLieAlgebraKey(type::Symbol, n::Int, F::String)
        @req type in [:gl, :so] "type must be either :gl or :so"
        @req n > 0 "n must be a positive integer"
        return new(type, n, F)
    end
end

struct DBLieAlgebraModuleKey
    modstring::String

    function DBLieAlgebraModuleKey(modstring::String)
        @req !isempty(modstring) "Module string must not be empty"
        return new(modstring)
    end
end

struct DBSmashProductKey
    L::DBLieAlgebraKey
    V::DBLieAlgebraModuleKey

    function DBSmashProductKey(Lk::DBLieAlgebraKey, Vk::DBLieAlgebraModuleKey)
        return new(Lk, Vk)
    end
end


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
        for (degs, save_basis, save_pbw) in [(deg:deg, true, true), (0:deg, false, true)]
            @vprintln :PBWDeformationsDatabase "Computing GlnGraphDeformBasis for degrees $(degs)..." # ProgressMeter needs a linebreak
            b = GlnGraphDeformBasis(sp, degs)
            @vprintln :PBWDeformationsDatabase " Done"
            if save_basis
                @vprint :PBWDeformationsDatabase "Saving GlnGraphDeformBasis..."
                save(joinpath(path, string_for_filename(b)) * file_ext, b; serializer=PBWDeformations.JSONSerializerNoRefs())
                @vprintln :PBWDeformationsDatabase " Done"
            end
            @vprint :PBWDeformationsDatabase "Computing PBW deformations for degrees $(degs)..."
            ms = all_pbwdeformations(sp, b)
            @vprintln :PBWDeformationsDatabase " Done"
            if save_pbw
                @vprint :PBWDeformationsDatabase "Saving PBW deformations..."
                save(joinpath(path, string_for_filename_pbwdeforms(b)) * file_ext, ms; serializer=PBWDeformations.JSONSerializerNoRefs())
                @vprintln :PBWDeformationsDatabase " Done"
            end
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

function load_glngraph_deform_bases(base_path::String, sp::SmashProductLie; degree_type::Symbol=:pure)
    @req degree_type in [:pure] "`degree_type` must be `:pure`"
    path = prepare_loading(base_path, sp)

    deg = 0
    bs = GlnGraphDeformBasis{elem_type(sp)}[]
    while (degs = (deg:deg); filepath = joinpath(path, string_for_filename(GlnGraphDeformBasis, sp, degs)) * file_ext; isfile(filepath))
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

function load_pbwdeformations(base_path::String, sp::SmashProductLie; degree_type::Symbol=:pure)
    @req degree_type in [:pure, :upto] "`degree_type` must be either `:pure` or `:upto`"
    path = prepare_loading(base_path, sp)

    deg = 0
    mss = Vector{DeformationMap{elem_type(sp)}}[]
    while (degs = degree_type == :pure ? (deg:deg) : (0:deg); filepath = joinpath(path, string_for_filename_pbwdeforms(sp, degs)) * file_ext; isfile(filepath))
        @vprint :PBWDeformationsDatabase "Found PBW deformations for degree $(degs). Loading..."
        push!(mss, Vector{DeformationMap{elem_type(sp)}}(load(filepath))::Vector{DeformationMap{elem_type(sp)}}) # see https://github.com/oscar-system/Oscar.jl/issues/3983
        @vprintln :PBWDeformationsDatabase " Done"
        deg += 1
    end
    return mss
end

function are_all_pbwdeformations_puredimensional(base_path::String, sp::SmashProductLie)
    pure_dims = length.(load_pbwdeformations(base_path, sp; degree_type=:pure))
    upto_dims = length.(load_pbwdeformations(base_path, sp; degree_type=:upto))
    return is_prefix_equal(cumsum(pure_dims), upto_dims)
end


################################################################################
#
# File path helpers
#
################################################################################

const file_ext = ".mrdi"

function DBLieAlgebraKey(type::Union{Symbol, String}, n::Int, F::Union{Field, Symbol, String})
    F_str = if F isa Field
        if F isa QQField
            "QQ"
        else
            error("currently only QQ is supported for the coefficient field")
        end
    elseif F isa Symbol
        string(F)
    elseif F isa String
        F
    end
    return DBLieAlgebraKey(Symbol(type), n, F_str)
end

function DBLieAlgebraKey(L::LinearLieAlgebra)
    type = type_short(L)
    n = L.n
    F = coefficient_ring(L)
    return DBLieAlgebraKey(type, n, F)
end

function string_for_filename(L::LinearLieAlgebra)
    return string_for_filename(DBLieAlgebraKey(L))
end

function string_for_filename(Lk::DBLieAlgebraKey)
    return "$(Lk.type)_$(Lk.n)_$(Lk.F)"
end

function type_short(L::LinearLieAlgebra)
    type = get_attribute(L, :type, nothing)
    if type == :general_linear
        return :gl
    elseif type == :special_orthogonal
        return :so
    end
    error("type_short not implemented for this type of Lie algebra")
end


function DBLieAlgebraModuleKey(V::LieAlgebraModule)
    return DBLieAlgebraModuleKey(module_as_string(V))
end

function string_for_filename(V::LieAlgebraModule)
    return string_for_filename(DBLieAlgebraModuleKey(V))
end

function string_for_filename(Vk::DBLieAlgebraModuleKey)
    return Vk.modstring
end

function module_as_string(V::LieAlgebraModule)
    if _is_standard_module(V)
        return "V"
    elseif ((fl, B) = _is_dual(V); fl)
        if _is_standard_module(B)
            return "DV"
        end
        return "D_" * module_as_string(B) * "_"
    elseif ((fl, Bs) = _is_direct_sum(V); fl)
        return "_" * join(module_as_string.(Bs), "_+_") * "_"
    elseif ((fl, Bs) = _is_tensor_product(V); fl)
        return "_" * join(module_as_string.(Bs), "_x_") * "_"
    elseif ((fl, B, k) = _is_exterior_power(V); fl)
        return "E$(k)_" * module_as_string(B) * "_"
    elseif ((fl, B, k) = _is_symmetric_power(V); fl)
        return "S$(k)_" * module_as_string(B) * "_"
    elseif ((fl, B, k) = _is_tensor_power(V); fl)
        return "T$(k)_" * module_as_string(B) * "_"
    end
    error("not implemented for this type of module")
end


function DBSmashProductKey(L::LinearLieAlgebra, V::LieAlgebraModule)
    return DBSmashProductKey(DBLieAlgebraKey(L), DBLieAlgebraModuleKey(V))
end

function DBSmashProductKey(sp::SmashProductLie)
    return DBSmashProductKey(base_lie_algebra(sp), base_module(sp))
end

function string_for_filename(sp::SmashProductLie)
    return string_for_filename(DBSmashProductKey(sp))
end

function string_for_filename(spk::DBSmashProductKey)
    return string_for_filename(spk.L) * "-" * string_for_filename(spk.V)
end


function string_for_filename(b::ArcDiagBasedDeformBasis)
    return string_for_filename(typeof(b), b.sp, b.degs)
end

function string_for_filename(T::Type{<:ArcDiagBasedDeformBasis}, sp::SmashProductLie, degs::AbstractVector{Int})
    return string_for_filename(T, DBSmashProductKey(sp), degs)
end

function string_for_filename(T::Type{<:ArcDiagDeformBasis}, spk::DBSmashProductKey, degs::AbstractVector{Int})
    return "ArcDiagDeformBasis-" * string_for_filename(spk) * "-" * join(degs, "_")
end

function string_for_filename(::Type{<:GlnGraphDeformBasis}, spk::DBSmashProductKey, degs::AbstractVector{Int})
    return "GlnGraphDeformBasis-" * string_for_filename(spk) * "-" * join(degs, "_")
end


function string_for_filename_pbwdeforms(b::ArcDiagBasedDeformBasis)
    return string_for_filename_pbwdeforms(b.sp, b.degs)
end

function string_for_filename_pbwdeforms(sp::SmashProductLie, degs::AbstractVector{Int})
    return string_for_filename_pbwdeforms(DBSmashProductKey(sp), degs)
end

function string_for_filename_pbwdeforms(spk::DBSmashProductKey, degs::AbstractVector{Int})
    return "PBWDeformations-" * string_for_filename(spk) * "-" * join(degs, "_")
end


function string_for_filename_setup(b::ArcDiagBasedDeformBasis)
    return string_for_filename_setup(b.sp)
end

function string_for_filename_setup(sp::SmashProductLie)
    return string_for_filename_setup(DBSmashProductKey(sp))
end

function string_for_filename_setup(spk::DBSmashProductKey)
    return "setup-" * string_for_filename(spk)
end


function string_for_path(b::ArcDiagBasedDeformBasis)
    return string_for_path(b.sp)
end

function string_for_path(sp::SmashProductLie)
    return string_for_path(DBSmashProductKey(sp))
end

function string_for_path(spk::DBSmashProductKey)
    return joinpath(string(spk.L.type), string_for_filename(spk.V))
end

end # module
