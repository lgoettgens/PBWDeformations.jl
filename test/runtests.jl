include("setup.jl")

# start tests
@test all(n -> isdefined(PBWDeformations, n), names(PBWDeformations))

include("ArcDiagram-test.jl")
include("DeformationBases-test.jl")
include("Pseudograph-test.jl")
include("SmashProductLie-test.jl")
include("SmashProductDeformLie-test.jl")
include("SmashProductPBWDeformLie-test.jl")
include("Util-test.jl")

DocMeta.setdocmeta!(
    PBWDeformations,
    :DocTestSetup,
    :(using PBWDeformations; using PBWDeformations.Oscar);
    recursive=true,
)
doctest(PBWDeformations)
