include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    include("Util-tets.jl")
    include("QuadraticAlgebra-test.jl")
    # include("PBWDeformations-test.jl")
    # include("AlgebraElement-test_old.jl")
    # include("QuadraticAlgebra-test.jl")
    # include("GroupAlgebra-test.jl")
    # include("SmashProductLie-test.jl")
    # include("SmashProductDeformLie-test.jl")
end
