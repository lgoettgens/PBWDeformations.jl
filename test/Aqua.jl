using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(
        PBWDeformations;
        ambiguities=false,          # recursive=false does not work
        persistent_tasks=false,     # does not work with Oscar version juggling
    )
end
