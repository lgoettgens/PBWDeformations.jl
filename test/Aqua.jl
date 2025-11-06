using Aqua

@testset "Aqua.jl" begin
    Aqua.test_all(
        PBWDeformations;
        ambiguities=false,          # recursive=false does not work
        persistent_tasks=false,     # does not work with Oscar version juggling
        piracies=(; treat_as_own=[Oscar.load_object, Oscar.Partition, Oscar.id_hom]) # TODO: remove when upstreamed
    )
end
