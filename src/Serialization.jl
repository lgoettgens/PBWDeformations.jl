function patch_oscar_serialization_namespace()
    Oscar.get_oscar_serialization_version() # call once to ensure the Oscar version is set

    push!(
        Oscar.oscar_serialization_version[],
        :PBWDeformations => ["https://github.com/lgoettgens/PBWDeformations.jl", Base.get_pkgversion_from_path(dirname(@__DIR__))],
    )
end
