include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    include("Util-tets.jl")
    include("QuadraticAlgebra-test.jl")
end
