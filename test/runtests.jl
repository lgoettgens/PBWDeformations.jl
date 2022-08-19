include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    include("Util-test.jl")
    include("ArcDiagram-test.jl")
    include("FreeAlgebra-test.jl")
    include("SmashProductLie-test.jl")
    include("SmashProductDeformLie-test.jl")

    DocMeta.setdocmeta!(PBWDeformations, :DocTestSetup, :(using PBWDeformations; using Oscar); recursive=true)
    doctest(PBWDeformations)
end
