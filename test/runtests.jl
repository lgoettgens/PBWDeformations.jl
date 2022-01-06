include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    include("Util-test.jl")
    include("FreeAlgebra-test.jl")
    include("SmashProductLie-test.jl")
    include("SmashProductDeformLie-test.jl")
end
