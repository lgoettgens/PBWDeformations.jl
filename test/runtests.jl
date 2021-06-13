include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    include("PBWDeformations-test.jl")
    include("SmashProductLie-test.jl")
    include("SmashProductSymmDeformLie-test.jl")
    include("AlgebraWithCommutators-test.jl")
end
