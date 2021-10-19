include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    include("PBWDeformations-test.jl")
    include("AlgebraElement-test.jl")
    include("QuadraticAlgebra-test.jl")
    include("GroupAlgebra-test.jl")
    include("SmashProductLie-test.jl")
    include("SmashProductDeformLie-test.jl")
end
