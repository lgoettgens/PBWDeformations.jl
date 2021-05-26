include("setup.jl")

@testset ExtendedTestSet "All PBWDeformations tests" begin
    @test 2 == PBWDeformations.add(1,1)
    # Write your tests here.
end
