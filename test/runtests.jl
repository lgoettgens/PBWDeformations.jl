include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    include("PBWDeformations-test.jl")
    include("QuadraticAlgebra-test.jl")
    include("SmashProductLie-test.jl")
    include("SmashProductSymmDeformLie-test.jl")
    include("SmashProductDeformLie-test.jl")
end
