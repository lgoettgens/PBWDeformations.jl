module TestExt

using PBWDeformations
using Test

using PBWDeformations.Oscar: Oscar, save, load

function PBWDeformations.test_save_load_roundtrip(func, path, original::T;
                                    params=nothing, check_func=nothing, kw...) where {T}
    # save and load from a file
    filename = joinpath(path, "original.json")
    save(filename, original; kw...)
    loaded = load(filename; params=params, kw...)

    @test loaded isa T
    func(loaded)

    # save and load from a file without saving references
    filename = joinpath(path, "original.json")
    save(filename, original; serializer=PBWDeformations.JSONSerializerNoRefs(), kw...)
    @test !any(line -> contains(line, "_refs"), eachline(filename))
    loaded = load(filename; params=params, serializer=PBWDeformations.JSONSerializerNoRefs(), kw...)

    @test loaded isa T
    func(loaded)

    # save and load from an IO buffer
    io = IOBuffer()
    save(io, original; kw...)
    seekstart(io)
    loaded = load(io; params=params, kw...)

    @test loaded isa T
    func(loaded)

    # save and load from an IO buffer, with prescribed type
    io = IOBuffer()
    save(io, original; kw...)
    seekstart(io)
    loaded = load(io; type=T, params=params, kw...)

    @test loaded isa T
    func(loaded)

    # test loading on a empty state
    save(filename, original; kw...)
    Oscar.reset_global_serializer_state()
    loaded = load(filename; params=params, kw...)
    @test loaded isa T

    isnothing(check_func) || @test check_func(loaded)
end

end # module
