include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    include("Util-test.jl")
    include("FreeAlgebra-test.jl")
    include("QuadraticAlgebra-test.jl")
end
