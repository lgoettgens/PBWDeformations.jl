include("setup.jl")

Oscar.with_unicode(true) do
    # short
    include("ArcDiagram-test.jl")
    include("DeformationBases-test.jl")
    include("LinearIndependence-test.jl")
    include("Pseudograph-test.jl")
    include("SmashProductLie-test.jl")
    include("SmashProductLieDeform-test.jl")
    include("SmashProductPBWDeformLie-test.jl")

    # long
    include("Aqua.jl")
    include("ModuleSimpleStructure-test.jl")

    if VERSION >= v"1.7-"
        DocMeta.setdocmeta!(
            PBWDeformations,
            :DocTestSetup,
            :(using PBWDeformations; using PBWDeformations.Oscar);
            recursive=true,
        )
        doctest(PBWDeformations)
    end
end
