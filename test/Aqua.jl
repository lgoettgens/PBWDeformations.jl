using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(
        PBWDeformations;
        ambiguities=false,          # recursive=false does not work
    )
end
